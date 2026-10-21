!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @brief This module provides a method of offloading I/O operations to dedicated writer PEs in order
!! to perform asynchronous I/O operations while the model PEs continue their computation.
!!
!! The writer PEs are split evenly between the domain tiles, so the number of writers must be a
!! multiple of the number of tiles (e.g. 6, 12, ... for a cubed sphere built by
!! mpp_define_cubic_domain). Writers (t-1)*n+1 to t*n write tile t, splitting it with a writer
!! layout.
!!
!! The writers can specify their own io layout, which determines whether they will write a single file
!! or if each writer will write out a portion of the data. While this allows for more flexibility,
!! using a 1,1 io layout with more than 1 PE will only slow down performance, since the data will be
!! gathered and written from the root writer PE.
!!
!!
!! The model side never waits on the writers during normal operation:
!! - Field data does not go through a model-side gather. Each model PE copies its compute domain
!!   into a buffer from a small ring and posts nonblocking sends of the parts each writer owns
!!   directly to those writers. A send only has to finish when its buffer is reused or the file
!!   is closed.
!! - Each operation (open, register, write, ...) is one packed control message from the model root
!!   to every writer, also sent without waiting from a ring of buffers.
!! Writers receive field data straight into place through MPI subarray types, and all writers
!! write at the same time. The number of simultaneously open files is set by max_files in
!! offloading_io_nml; offloading_io_end ends the writer session.
module offloading_io_mod
#ifdef use_libMPI
!! TODO mpi calls should all go through mpp and have a non f08 interface ideally
  use mpi_f08, only: MPI_Comm, MPI_Request, MPI_Datatype, MPI_Isend, MPI_Irecv, MPI_Recv, MPI_Send, MPI_Wait, &
                     MPI_Waitall, MPI_Testall, MPI_Type_create_subarray, MPI_Type_commit, MPI_Type_free, &
                     MPI_INTEGER, MPI_REAL4, MPI_REAL8, MPI_ORDER_FORTRAN, MPI_STATUS_IGNORE, MPI_DATATYPE_NULL, &
                     MPI_STATUSES_IGNORE, MPI_REQUEST_NULL
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_get_current_pelist, mpp_declare_pelist, mpp_set_current_pelist, &
                    mpp_sync, mpp_sum, mpp_error, FATAL, input_nml_file
  use fms_mod, only: check_nml_error
  use mpp_domains_mod
  use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, register_axis, register_field, write_data, &
                         register_global_attribute
  use platform_mod, only: r4_kind, r8_kind, i4_kind, i8_kind

  implicit none
  private

  integer :: max_files = 10
  namelist /offloading_io_nml/ max_files

  !> opcode enumerations to send the writers
  integer, parameter :: IO_OP_OPEN = 1, IO_OP_REGISTER_AXIS = 2, IO_OP_REGISTER_FIELD = 3
  integer, parameter :: IO_OP_WRITE = 4, IO_OP_CLOSE = 5, IO_OP_GLOBAL_ATTRIBUTE = 6, IO_OP_SHUTDOWN = 7
  integer, parameter :: ATTRIBUTE_R4 = 1, ATTRIBUTE_R8 = 2, ATTRIBUTE_I4 = 3, ATTRIBUTE_I8 = 4, ATTRIBUTE_STR = 5
  integer, parameter :: FIELD_R4 = 1, FIELD_R8 = 2
  !> These are arbitrary values used to tag the data flow in mpi communication
  !! protocol -> sending the opcodes, data -> sending the actual data to write
  integer, parameter :: protocol_tag = 24001, data_tag = 24002, string_length = 255
  integer, parameter :: default_max_pending_writes = 4
  !> Strings are packed into integers, 4 characters per word for 4 byte integers.
  integer, parameter :: string_words = (string_length + storage_size(0)/8) / (storage_size(0)/8)
  !> A control message holds an operation and at most two strings plus a few integer words.
  integer, parameter :: control_length = 2*string_words + 8
  !> Number of control messages the model root can have in flight to each writer.
  integer, parameter :: control_slots = 64
  !> Columns of the model layout table: tile, isc, iec, jsc, jec, tile x size, tile y size,
  !! number of tiles, writer layout (2), writer io_layout (2)
  integer, parameter :: layout_info_size = 12
  !> Columns of the writer layout table: tile, isc, iec, jsc, jec
  integer, parameter :: writer_info_size = 5
  integer, allocatable :: all_pes(:), model_pes(:), writer_pes(:)
  integer, allocatable :: writer_ranks(:)
  integer :: model_root_rank = -1
  integer :: global_size(2) = 0 !< size of each tile
  integer :: ntiles = 0
  integer :: writer_layout(2) = 0, writer_io_layout(2) = 0 !< decomposition of each tile among its writers
  logical :: metadata_ready = .false., session_finished = .false.
  logical, allocatable :: file_is_open(:)
  type(MPI_Comm) :: offload_comm

  ! Control message being packed (model root) or unpacked (writers).
  integer :: message(control_length), message_position = 0

  ! Model root state: ring of control messages in flight to the writers.
  integer :: next_control = 1
  integer, allocatable, asynchronous :: control_buffers(:,:)
  type(MPI_Request), allocatable, asynchronous :: control_requests(:,:)

  ! Model PE state: compute bounds, the writers that own parts of it and which parts,
  ! and the ring of field data in flight (send_requests(writer, slot)).
  integer :: compute_bounds(4) = 0, next_send = 1
  integer, allocatable :: dest_ranks(:)
  type(MPI_Datatype), allocatable :: send_types(:,:)
  real(kind=r4_kind), allocatable, asynchronous :: send_buffers_r4(:,:,:)
  real(kind=r8_kind), allocatable, asynchronous :: send_buffers_r8(:,:,:)
  type(MPI_Request), allocatable, asynchronous :: send_requests(:,:)

  ! Writer PE state: its part of the tile, the model ranks that overlap it and where their parts go.
  type(domain2d), save :: writer_domain
  type(FmsNetcdfDomainFile_t), allocatable :: writer_fileobjs(:)
  integer, allocatable :: writer_naxes(:), writer_axis_directions(:,:)
  character(len=string_length), allocatable :: writer_axis_names(:,:)
  integer, allocatable :: source_ranks(:)
  type(MPI_Datatype), allocatable :: recv_types(:,:)
  type(MPI_Request), allocatable, asynchronous :: recv_requests(:)
  real(kind=r4_kind), allocatable, asynchronous :: block_data_r4(:,:)
  real(kind=r8_kind), allocatable, asynchronous :: block_data_r8(:,:)

  public :: setup_offloading_pelists, send_domain_metadata, recv_domain_metadata, run_offload_writer
  public :: open_offloaded_file, close_offloaded_file, offloading_io_end
  public :: register_axis_offloaded_file, register_axes_offloaded_file
  public :: register_field_offloaded_file, write_field_offloaded_file, global_metadata_offload

  !> Interface for writing fields to an offloaded file, supporting both real(kind=r4_kind) and real(kind=r8_kind) data
  !! TODO: more data types
  interface write_field_offloaded_file
    module procedure write_field_offloaded_file_r4
    module procedure write_field_offloaded_file_r8
  end interface

contains

  !> Collective on the initial pelist; reserve its last nwriters PEs for the writers.
  !! nwriters must be a multiple of the number of tiles in the model domain (6 for a cubed sphere).
  subroutine setup_offloading_pelists(global_pelist, model_pelist, writer_pelist, nwriters)
    integer, allocatable, intent(out) :: global_pelist(:), model_pelist(:), writer_pelist(:)
    integer, intent(in), optional :: nwriters !< number of writer PEs (default 1)
    integer :: n, comm_id, w, io_status

    if (allocated(all_pes)) call mpp_error(FATAL, "offloading: pelists already initialized")
    n = 1
    if (present(nwriters)) n = nwriters
    if (n < 1) call mpp_error(FATAL, "offloading: at least one writer PE is required")
    if (mpp_npes() <= n) call mpp_error(FATAL, "offloading: at least one model PE is required besides the writers")
    read(input_nml_file, nml=offloading_io_nml, iostat=io_status)
    io_status = check_nml_error(io_status, "offloading_io_nml")
    if (max_files < 1) call mpp_error(FATAL, "offloading: max_files must be positive")
    allocate(file_is_open(max_files))
    file_is_open = .false.
    allocate(all_pes(mpp_npes()))
    call mpp_get_current_pelist(all_pes)
    model_pes = all_pes(:size(all_pes)-n)
    writer_pes = all_pes(size(all_pes)-n+1:)
    call mpp_declare_pelist(all_pes)
    call mpp_declare_pelist(model_pes)
    call mpp_declare_pelist(writer_pes)
    call mpp_set_current_pelist(all_pes, no_sync=.true.)
    call mpp_get_current_pelist(all_pes, commID=comm_id)
    offload_comm%mpi_val = comm_id
    model_root_rank = rank_of(model_pes(1))
    allocate(writer_ranks(n))
    do w = 1, n
      writer_ranks(w) = rank_of(writer_pes(w))
    enddo
    global_pelist = all_pes
    model_pelist = model_pes
    writer_pelist = writer_pes
  end subroutine

  !> All model PEs call this once, paired with recv_domain_metadata on every writer.
  !! Each model PE must own exactly one tile, and the number of writer PEs must be a multiple of the
  !! number of tiles. All optional arguments must be the same on every model PE.
  subroutine send_domain_metadata(model_domain, max_pending_writes, writer_layout, writer_io_layout)
    type(domain2d), intent(inout) :: model_domain
    integer, intent(in), optional :: max_pending_writes !< field writes each model PE can have in
                                                        !! flight (default 4)
    integer, intent(in), optional :: writer_layout(2) !< layout of each tile's writers
                                                      !! (default from mpp_define_layout)
    integer, intent(in), optional :: writer_io_layout(2) !< io_layout of each tile's writers (default
                                                         !! (1,1), one file per tile); must divide
                                                         !! writer_layout
    integer :: info(layout_info_size), tile_id(1), nslots, ntile_count, writers_per_tile

    call require_model()
    if (metadata_ready) call mpp_error(FATAL, "offloading: metadata already sent")
    ntile_count = mpp_get_ntile_count(model_domain)
    if (mod(size(writer_pes), ntile_count) /= 0) &
      call mpp_error(FATAL, "offloading: the number of writer PEs must be a multiple of the number of domain tiles")
    if (mpp_get_current_ntile(model_domain) /= 1) &
      call mpp_error(FATAL, "offloading: each model PE must own exactly one tile")
    writers_per_tile = size(writer_pes) / ntile_count
    nslots = default_max_pending_writes
    if (present(max_pending_writes)) nslots = max_pending_writes
    if (nslots < 1) call mpp_error(FATAL, "offloading: max_pending_writes must be positive")
    tile_id = mpp_get_tile_id(model_domain)
    info(1) = tile_id(1)
    call mpp_get_compute_domain(model_domain, info(2), info(3), info(4), info(5))
    call mpp_get_global_domain(model_domain, xsize=info(6), ysize=info(7))
    info(8) = ntile_count
    if (present(writer_layout)) then
      info(9:10) = writer_layout
    else
      call mpp_define_layout((/1, info(6), 1, info(7)/), writers_per_tile, info(9:10))
    endif
    info(11:12) = 1
    if (present(writer_io_layout)) info(11:12) = writer_io_layout
    call mpp_set_current_pelist(all_pes, no_sync=.true.)
    call mpp_broadcast_domain(model_domain)
    call exchange_layout(info)
    allocate(send_buffers_r4(info(3)-info(2)+1, info(5)-info(4)+1, nslots))
    allocate(send_buffers_r8(info(3)-info(2)+1, info(5)-info(4)+1, nslots))
    allocate(send_requests(size(dest_ranks), nslots))
    send_requests = MPI_REQUEST_NULL
    next_send = 1
    if (mpp_pe() == model_pes(1)) then
      allocate(control_buffers(control_length, control_slots))
      allocate(control_requests(size(writer_pes), control_slots))
      control_requests = MPI_REQUEST_NULL
      next_control = 1
    endif
    call mpp_set_current_pelist(model_pes, no_sync=.true.)
    metadata_ready = .true.
  end subroutine

  !> Called by every writer PE. global_domain_size returns the size of each tile.
  subroutine recv_domain_metadata(model_domain, global_domain_size)
    type(domain2d), intent(inout) :: model_domain
    integer, intent(out) :: global_domain_size(2)
    integer :: info(layout_info_size)

    call require_setup()
    if (.not. any(writer_pes == mpp_pe())) call mpp_error(FATAL, "offloading: receiver must be a writer PE")
    if (metadata_ready) call mpp_error(FATAL, "offloading: metadata already received")
    call mpp_set_current_pelist(all_pes, no_sync=.true.)
    call mpp_define_null_domain(model_domain)
    call mpp_broadcast_domain(model_domain)
    info = 0
    call exchange_layout(info)
    allocate(writer_fileobjs(max_files), writer_naxes(max_files), &
         writer_axis_directions(2,max_files), writer_axis_names(2,max_files))
    writer_naxes = 0
    writer_axis_directions = 0
    writer_axis_names = ""
    global_domain_size = global_size
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    metadata_ready = .true.
  end subroutine

  !> Collective on all PEs. Shares the model PEs' tiles and compute bounds, defines the writers'
  !! decomposition of the tiles, and sets up the MPI types that move each overlap between a model
  !! PE and a writer.
  subroutine exchange_layout(info)
    integer, intent(in) :: info(layout_info_size) !< this PE's layout information (zeros on writers)
    integer, allocatable :: table(:,:), writer_table(:,:), global_indices(:,:), layouts(:,:)
    integer, allocatable :: pe_start(:), pe_end(:)
    integer :: n, my_model, my_writer, ncells, writers_per_tile, overlap(4), no_contact(0)
    integer :: overlap_cells

    allocate(table(layout_info_size, size(model_pes)))
    table = 0
    my_model = findloc(model_pes, mpp_pe(), 1)
    if (my_model > 0) table(:, my_model) = info
    call mpp_sum(table, size(table))
    global_size = table(6:7, 1)
    ntiles = table(8, 1)
    writer_layout = table(9:10, 1)
    writer_io_layout = table(11:12, 1)
    do n = 1, size(model_pes)
      if (any(table(6:12,n) /= table(6:12,1))) &
        call mpp_error(FATAL, "offloading: tile size, writer layout and io_layout must agree on all model PEs")
      if (table(1,n) < 1 .or. table(1,n) > ntiles) call mpp_error(FATAL, "offloading: invalid model tile id")
    enddo
    do n = 1, ntiles
      if (.not. any(table(1,:) == n)) call mpp_error(FATAL, "offloading: a tile has no model PEs")
    enddo
    writers_per_tile = size(writer_pes) / ntiles
    if (writers_per_tile * ntiles /= size(writer_pes)) &
      call mpp_error(FATAL, "offloading: the number of writer PEs must be a multiple of the number of domain tiles")
    if (any(writer_layout < 1) .or. product(writer_layout) /= writers_per_tile) &
      call mpp_error(FATAL, "offloading: writer_layout must use all the writers of a tile")
    if (any(writer_io_layout < 1) .or. any(mod(writer_layout, max(writer_io_layout, 1)) /= 0)) &
      call mpp_error(FATAL, "offloading: writer_io_layout must divide writer_layout")
    if (any(writer_layout > global_size)) call mpp_error(FATAL, "offloading: writer_layout is larger than a tile")

    ! Writers define their decomposition: tile t goes to writers (t-1)*writers_per_tile+1 onward.
    allocate(writer_table(writer_info_size, size(writer_pes)))
    writer_table = 0
    my_writer = findloc(writer_pes, mpp_pe(), 1)
    if (my_writer > 0) then
      allocate(global_indices(4, ntiles), layouts(2, ntiles), pe_start(ntiles), pe_end(ntiles))
      do n = 1, ntiles
        global_indices(:,n) = (/1, global_size(1), 1, global_size(2)/)
        layouts(:,n) = writer_layout
        pe_start(n) = writer_pes((n-1)*writers_per_tile + 1)
        pe_end(n) = writer_pes(n*writers_per_tile)
      enddo
      call mpp_set_current_pelist(writer_pes, no_sync=.true.)
      ! No contacts: the domain only describes how the writers split and name the tile files.
      call mpp_define_mosaic(global_indices, layouts, writer_domain, ntiles, 0, no_contact, no_contact, &
                             no_contact, no_contact, no_contact, no_contact, no_contact, no_contact, no_contact, &
                             no_contact, pe_start, pe_end, name="offload writer mosaic")
      call mpp_define_io_domain(writer_domain, writer_io_layout)
      writer_table(1, my_writer) = (my_writer - 1) / writers_per_tile + 1
      call mpp_get_compute_domain(writer_domain, writer_table(2, my_writer), writer_table(3, my_writer), &
                                  writer_table(4, my_writer), writer_table(5, my_writer))
      call mpp_set_current_pelist(all_pes, no_sync=.true.)
    endif
    call mpp_sum(writer_table, size(writer_table))

    if (my_model > 0) then
      ! Model PE: one send type per writer that owns part of this compute domain.
      call find_overlaps(table(:, my_model), writer_table, writer_pes, dest_ranks, send_types)
      ncells = (info(3)-info(2)+1) * (info(5)-info(4)+1)
      compute_bounds = info(2:5)
    else
      ! Writer: one receive type per model PE that owns part of this writer's block.
      call find_overlaps(writer_table(:, my_writer), table, model_pes, source_ranks, recv_types)
      allocate(recv_requests(size(source_ranks)))
      recv_requests = MPI_REQUEST_NULL
      allocate(block_data_r4(writer_table(3,my_writer)-writer_table(2,my_writer)+1, &
                 writer_table(5,my_writer)-writer_table(4,my_writer)+1))
      allocate(block_data_r8(writer_table(3,my_writer)-writer_table(2,my_writer)+1, &
                 writer_table(5,my_writer)-writer_table(4,my_writer)+1))
      ncells = size(block_data_r4)
    endif
    if (overlap_cells /= ncells) call mpp_error(FATAL, "offloading: model and writer domains do not match")

  contains

    !> Find the PEs in others (same tile) whose bounds overlap, and build a subarray type for each
    !! overlap, relative to my own bounds.
    subroutine find_overlaps(mine, others, other_pes, ranks, types)
      integer, intent(in) :: mine(:), others(:,:), other_pes(:)
      integer, allocatable, intent(out) :: ranks(:)
      type(MPI_Datatype), allocatable, intent(out) :: types(:,:)
      integer :: m, k

      allocate(ranks(0), types(2,0))
      overlap_cells = 0
      do m = 1, size(other_pes)
        if (others(1,m) /= mine(1)) cycle
        overlap = (/max(mine(2), others(2,m)), min(mine(3), others(3,m)), &
                    max(mine(4), others(4,m)), min(mine(5), others(5,m))/)
        if (overlap(1) > overlap(2) .or. overlap(3) > overlap(4)) cycle
        k = size(ranks) + 1
        ranks = (/ranks, rank_of(other_pes(m))/)
        types = reshape([types, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL], [2, k])
        call MPI_Type_create_subarray(2, (/mine(3)-mine(2)+1, mine(5)-mine(4)+1/), &
                                      (/overlap(2)-overlap(1)+1, overlap(4)-overlap(3)+1/), &
                                      (/overlap(1)-mine(2), overlap(3)-mine(4)/), MPI_ORDER_FORTRAN, &
                                      MPI_REAL4, types(1,k))
        call MPI_Type_commit(types(1,k))
        call MPI_Type_create_subarray(2, (/mine(3)-mine(2)+1, mine(5)-mine(4)+1/), &
                    (/overlap(2)-overlap(1)+1, overlap(4)-overlap(3)+1/), &
                    (/overlap(1)-mine(2), overlap(3)-mine(4)/), MPI_ORDER_FORTRAN, &
                    MPI_REAL8, types(2,k))
        call MPI_Type_commit(types(2,k))
        overlap_cells = overlap_cells + (overlap(2)-overlap(1)+1) * (overlap(4)-overlap(3)+1)
      enddo
    end subroutine
  end subroutine

  subroutine open_offloaded_file(filename, fileobj)
    character(len=*), intent(in) :: filename
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: file_id

    call require_model()
    if (.not. metadata_ready) call mpp_error(FATAL, "offloading: send domain metadata before opening a file")
    if (session_finished) call mpp_error(FATAL, "offloading: session already finished")
    call check_string(filename)
    file_id = find_free_file_id()
    if (file_id == 0) call mpp_error(FATAL, "offloading: max_files exceeded")
    file_is_open(file_id) = .true.
    fileobj%offloading_obj_in%id = file_id
    call start_message(IO_OP_OPEN)
    call pack_words((/file_id/))
    call pack_string(filename)
    call post_message()
  end subroutine

  !> axis_size is the size of one tile in the cart_name direction.
  subroutine register_axis_offloaded_file(axis_name, cart_name, axis_size, fileobj)
    character(len=*), intent(in) :: axis_name, cart_name
    integer, intent(in) :: axis_size
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: direction, file_id

    file_id = require_file(fileobj)
    call check_string(axis_name)
    select case (cart_name)
    case ("x")
      direction = 1
    case ("y")
      direction = 2
    case default
      call mpp_error(FATAL, "offloading: axis cart_name must be x or y")
    end select
    if (axis_size /= global_size(direction)) call mpp_error(FATAL, "offloading: axis size differs from domain")
    call start_message(IO_OP_REGISTER_AXIS)
    call pack_words((/file_id/))
    call pack_string(axis_name)
    call pack_string(cart_name)
    call pack_words((/axis_size/))
    call post_message()
  end subroutine

  subroutine register_axes_offloaded_file(axis_names, cart_names, axis_sizes, domain, fileobj)
    character(len=*), intent(in) :: axis_names(:), cart_names(:)
    integer, intent(in) :: axis_sizes(:)
    type(domain2d), intent(in) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: i, nx, ny, file_id

    file_id = require_file(fileobj)
    if (size(axis_names) /= size(cart_names) .or. size(axis_names) /= size(axis_sizes)) &
      call mpp_error(FATAL, "offloading: axis arrays must have equal sizes")
    call mpp_get_global_domain(domain, xsize=nx, ysize=ny)
    if (any(global_size /= (/nx, ny/))) call mpp_error(FATAL, "offloading: axis domain differs from metadata")
    do i = 1, size(axis_names)
      call register_axis_offloaded_file(axis_names(i), cart_names(i), axis_sizes(i), fileobj)
    enddo
  end subroutine

  !> Dimension indices refer to registration order, and must describe x then y.
  !! TODO, should field type/kind argument be required? r8 default currently if not given
  subroutine register_field_offloaded_file(field_name, field_dims, field_type, fileobj)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: field_dims(:)
    character(len=*), intent(in), optional :: field_type
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: field_kind
    integer :: file_id

    file_id = require_file(fileobj)
    call check_string(field_name)
    if (size(field_dims) /= 2) call mpp_error(FATAL, "offloading: only 2-D fields are supported")
    if (any(field_dims < 1) .or. any(field_dims > 2)) &
      call mpp_error(FATAL, "offloading: field dimension indices must be 1 or 2")
    field_kind = FIELD_R8
    if (present(field_type)) then
      select case (trim(field_type))
      case ("float")
        field_kind = FIELD_R4
      case ("double")
        field_kind = FIELD_R8
      case default
        call mpp_error(FATAL, "offloading: field_type must be float or double")
      end select
    endif
    call start_message(IO_OP_REGISTER_FIELD)
    call pack_words((/file_id/))
    call pack_string(field_name)
    call pack_words(field_dims)
    call pack_words((/field_kind/))
    call post_message()
  end subroutine

  !> Register a scalar global attribute (r4, r8, i4, i8 or string) on the offloaded file.
  !! Only the model root's value is sent, and every tile file gets it.
  !! Can be called at any time while the file is open.
  subroutine global_metadata_offload(attribute_name, attribute_value, fileobj)
    character(len=*), intent(in) :: attribute_name
    class(*), intent(in) :: attribute_value
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: file_id

    file_id = require_file(fileobj)
    call check_string(attribute_name)
    call start_message(IO_OP_GLOBAL_ATTRIBUTE)
    call pack_words((/file_id/))
    call pack_string(attribute_name)
    select type (attribute_value)
    type is (real(kind=r4_kind))
      call pack_words((/ATTRIBUTE_R4/))
      call pack_scalar(attribute_value)
    type is (real(kind=r8_kind))
      call pack_words((/ATTRIBUTE_R8/))
      call pack_scalar(attribute_value)
    type is (integer(kind=i4_kind))
      call pack_words((/ATTRIBUTE_I4/))
      call pack_scalar(attribute_value)
    type is (integer(kind=i8_kind))
      call pack_words((/ATTRIBUTE_I8/))
      call pack_scalar(attribute_value)
    type is (character(len=*))
      if (len_trim(attribute_value) > string_length) &
        call mpp_error(FATAL, "offloading: string attribute values must contain at most 255 characters")
      call pack_words((/ATTRIBUTE_STR, len_trim(attribute_value)/))
      call pack_string(attribute_value)
    class default
      call mpp_error(FATAL, "offloading: unsupported global attribute type for "//trim(attribute_name))
    end select
    call post_message()
  end subroutine

  !> Collective on model PEs, but no PE waits on the others or on the writer. Field data can cover
  !! the compute or the data domain; halos are not sent. Returns once the data is copied into a
  !! send buffer, so the caller can reuse field_data right away.
  subroutine write_field_offloaded_file_r4(field_name, field_data, domain, fileobj)
    character(len=*), intent(in) :: field_name
    real(kind=r4_kind), intent(in) :: field_data(:,:)
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: ioff, joff, nx, ny, k, d
    integer :: file_id

    file_id = require_file(fileobj)
    call check_string(field_name)
    !! TODO, this might only be needed when writer PEs > 1. can probably optimize for single writer case
    call get_field_offsets(shape(field_data), domain, ioff, joff, nx, ny)

    ! Reuse the oldest buffer; this only waits if all max_pending_writes sends are still in flight.
    k = next_send
    call MPI_Waitall(size(dest_ranks), send_requests(:,k), MPI_STATUSES_IGNORE)
    next_send = mod(next_send, size(send_requests, 2)) + 1
    send_buffers_r4(:,:,k) = field_data(ioff+1:ioff+nx, joff+1:joff+ny)
    do d = 1, size(dest_ranks)
      call MPI_Isend(send_buffers_r4(:,:,k), 1, send_types(1,d), dest_ranks(d), data_tag, offload_comm, &
                     send_requests(d,k))
    enddo
    call start_message(IO_OP_WRITE)
    call pack_words((/file_id/))
    call pack_string(field_name)
    call pack_words((/FIELD_R4/))
    call post_message()
  end subroutine

  subroutine write_field_offloaded_file_r8(field_name, field_data, domain, fileobj)
    character(len=*), intent(in) :: field_name
    real(kind=r8_kind), intent(in) :: field_data(:,:)
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: ioff, joff, nx, ny, k, d
    integer :: file_id

    file_id = require_file(fileobj)
    call check_string(field_name)
    call get_field_offsets(shape(field_data), domain, ioff, joff, nx, ny)

    k = next_send
    call MPI_Waitall(size(dest_ranks), send_requests(:,k), MPI_STATUSES_IGNORE)
    next_send = mod(next_send, size(send_requests, 2)) + 1
    send_buffers_r8(:,:,k) = field_data(ioff+1:ioff+nx, joff+1:joff+ny)
    do d = 1, size(dest_ranks)
      call MPI_Isend(send_buffers_r8(:,:,k), 1, send_types(2,d), dest_ranks(d), data_tag, offload_comm, &
                     send_requests(d,k))
    enddo
    call start_message(IO_OP_WRITE)
    call pack_words((/file_id/))
    call pack_string(field_name)
    call pack_words((/FIELD_R8/))
    call post_message()
  end subroutine

  !> Figures out where to write the field data within the domain
  subroutine get_field_offsets(field_shape, domain, ioff, joff, nx, ny)
    integer, intent(in) :: field_shape(2)
    type(domain2d), intent(inout) :: domain
    integer, intent(out) :: ioff, joff
    integer, intent(out) :: nx, ny
    integer :: global_nx, global_ny, isc, iec, jsc, jec, isd, ied, jsd, jed

    call mpp_get_global_domain(domain, xsize=global_nx, ysize=global_ny)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)
    if (any(global_size /= (/global_nx, global_ny/)) .or. any(compute_bounds /= (/isc, iec, jsc, jec/))) &
      call mpp_error(FATAL, "offloading: field domain differs from metadata")
    nx = iec - isc + 1
    ny = jec - jsc + 1
    if (all(field_shape == (/iec-isc+1, jec-jsc+1/))) then
      ioff = 0
      joff = 0
    else if (all(field_shape == (/ied-isd+1, jed-jsd+1/))) then
      ioff = isc - isd
      joff = jsc - jsd
    else
      call mpp_error(FATAL, "offloading: field data must cover the compute or data domain")
    endif
  end subroutine

  !> Close one file handle, or close the legacy file and stop the writer session.
  subroutine close_offloaded_file(fileobj)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer :: file_id

    file_id = require_file(fileobj)
    call close_file_slot(file_id)
    fileobj%offloading_obj_in%id = 0
  end subroutine

  subroutine close_file_slot(file_id)
    integer, intent(in) :: file_id
    integer :: w

    if (file_id < 1 .or. file_id > max_files) call mpp_error(FATAL, "offloading: invalid file id")
    if (.not. file_is_open(file_id)) call mpp_error(FATAL, "offloading: no open file")
    do w = 1, size(send_requests, 2)
      call MPI_Waitall(size(dest_ranks), send_requests(:,w), MPI_STATUSES_IGNORE)
    enddo
    call start_message(IO_OP_CLOSE)
    call pack_words((/file_id/))
    call post_message()
    file_is_open(file_id) = .false.
  end subroutine

  subroutine offloading_io_end()
    integer :: acknowledgement, w, d, file_id

    call require_model()
    if (session_finished) call mpp_error(FATAL, "offloading: session already finished")
    do file_id = 1, max_files
      if (file_is_open(file_id)) call close_file_slot(file_id)
    enddo
    call start_message(IO_OP_SHUTDOWN)
    call post_message()
    if (mpp_pe() == model_pes(1)) then
      do w = 1, control_slots
        call MPI_Waitall(size(writer_pes), control_requests(:,w), MPI_STATUSES_IGNORE)
      enddo
      do w = 1, size(writer_pes)
        call MPI_Recv(acknowledgement, 1, MPI_INTEGER, writer_ranks(w), protocol_tag, offload_comm, &
                      MPI_STATUS_IGNORE)
        if (acknowledgement /= IO_OP_SHUTDOWN) call mpp_error(FATAL, "offloading: invalid shutdown acknowledgement")
      enddo
      deallocate(control_buffers, control_requests)
    endif
    call mpp_sync()
    do w = 1, size(send_types, 2)
      do d = 1, 2
        call MPI_Type_free(send_types(d,w))
      enddo
    enddo
    deallocate(send_buffers_r4, send_buffers_r8, send_requests, dest_ranks, send_types)
    deallocate(file_is_open)
    session_finished = .true.
  end subroutine

  !> Enters the main write loop for offloaded writer PE's to receive data and control messages from the model PEs
  !! To be called by every writer after recv_domain_metadata.
  subroutine run_offload_writer()
    integer :: operation, file_id, file_id_words(1), axis_size(1), direction, field_dims(2), field_kind(1), n
    character(len=string_length) :: name, cart

    ! checks for proper setup (setup_offloading_pelists,send/recv_domain_metadata)
    call require_setup()
    ! error checking for module state
    if (.not. any(writer_pes == mpp_pe())) call mpp_error(FATAL, "offloading: consumer must run on a writer PE")
    if (.not. metadata_ready) call mpp_error(FATAL, "offloading: receive metadata before running writer")
    if (session_finished) call mpp_error(FATAL, "offloading: writer session already finished")
    ! switch to the writers pelist if not already
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    ! main loop for io writers, receive operations from model PEs and handle them accordingly
    ! operations are "packed" into a single message for less overhead
    do
      operation = recv_message()
      select case (operation)
      case (IO_OP_OPEN)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call check_file_id(file_id)
        if (file_is_open(file_id)) call mpp_error(FATAL, "offloading: file already open")
        call unpack_string(name)
        if (.not. open_file(writer_fileobjs(file_id), trim(name), "overwrite", domain=writer_domain)) &
          call mpp_error(FATAL, "offloading: unable to open "//trim(name))
        file_is_open(file_id) = .true.
        writer_naxes(file_id) = 0
        writer_axis_directions(:,file_id) = 0
      case (IO_OP_REGISTER_AXIS)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call unpack_string(name)
        call unpack_string(cart)
        call unpack_words(axis_size)
        if (writer_naxes(file_id) == 2) call mpp_error(FATAL, "offloading: only two axes are supported")
        select case (trim(cart))
        case ("x")
          direction = 1
        case ("y")
          direction = 2
        case default
          call mpp_error(FATAL, "offloading: invalid axis direction received")
        end select
        if (axis_size(1) /= global_size(direction)) call mpp_error(FATAL, "offloading: invalid axis size received")
        if (any(writer_axis_directions(:writer_naxes(file_id),file_id) == direction)) &
          call mpp_error(FATAL, "offloading: duplicate axis direction")
        writer_naxes(file_id) = writer_naxes(file_id) + 1
        writer_axis_names(writer_naxes(file_id),file_id) = name
        writer_axis_directions(writer_naxes(file_id),file_id) = direction
        call register_axis(writer_fileobjs(file_id), trim(name), trim(cart))
      case (IO_OP_REGISTER_FIELD)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call unpack_string(name)
        call unpack_words(field_dims)
        call unpack_words(field_kind)
        if (any(field_dims < 1) .or. any(field_dims > writer_naxes(file_id))) &
          call mpp_error(FATAL, "offloading: field refers to unregistered axes")
        if (any(writer_axis_directions(field_dims,file_id) /= (/1,2/))) &
          call mpp_error(FATAL, "offloading: field dimensions must describe x then y")
        if (field_kind(1) == FIELD_R4) then
          call register_field(writer_fileobjs(file_id), trim(name), "float", writer_axis_names(field_dims,file_id))
        else if (field_kind(1) == FIELD_R8) then
          call register_field(writer_fileobjs(file_id), trim(name), "double", writer_axis_names(field_dims,file_id))
        else
          call mpp_error(FATAL, "offloading: invalid field type received")
        endif
      case (IO_OP_WRITE)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call unpack_string(name)
        call unpack_words(field_kind)
        select case (field_kind(1))
        case (FIELD_R4)
          !! TODO: can probably send the whole array at once instead of looping over source ranks
          do n = 1, size(source_ranks)
            call MPI_Irecv(block_data_r4, 1, recv_types(1,n), source_ranks(n), data_tag, offload_comm, &
                           recv_requests(n))
          enddo
          call MPI_Waitall(size(recv_requests), recv_requests, MPI_STATUSES_IGNORE)
          call write_data(writer_fileobjs(file_id), trim(name), block_data_r4)
        case (FIELD_R8)
          do n = 1, size(source_ranks)
            call MPI_Irecv(block_data_r8, 1, recv_types(2,n), source_ranks(n), data_tag, offload_comm, &
                           recv_requests(n))
          enddo
          call MPI_Waitall(size(recv_requests), recv_requests, MPI_STATUSES_IGNORE)
          call write_data(writer_fileobjs(file_id), trim(name), block_data_r8)
        case default
          call mpp_error(FATAL, "offloading: invalid field type received")
        end select
      case (IO_OP_GLOBAL_ATTRIBUTE)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call recv_global_attribute(writer_fileobjs(file_id))
      case (IO_OP_CLOSE)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call close_file(writer_fileobjs(file_id))
        file_is_open(file_id) = .false.
        writer_naxes(file_id) = 0
      case (IO_OP_SHUTDOWN)
        do file_id = 1, max_files
          if (file_is_open(file_id)) call close_file(writer_fileobjs(file_id))
        enddo
        session_finished = .true.
        call MPI_Send(IO_OP_SHUTDOWN, 1, MPI_INTEGER, model_root_rank, protocol_tag, offload_comm)
        exit
      case default
        call mpp_error(FATAL, "offloading: unknown operation received")
      end select
    enddo
    call mpp_deallocate_domain(writer_domain)
    do n = 1, size(recv_types, 2)
      call MPI_Type_free(recv_types(1,n))
      call MPI_Type_free(recv_types(2,n))
    enddo
    deallocate(source_ranks, recv_types, recv_requests, block_data_r4, block_data_r8)
    deallocate(writer_fileobjs, writer_naxes, writer_axis_directions, writer_axis_names)
    deallocate(file_is_open)
  end subroutine

  !> Writer half of global_metadata_offload; unpacks the rest of the current message.
  subroutine recv_global_attribute(fileobj)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    character(len=string_length) :: name, str_value
    integer :: value_type(1), str_len(1)
    real(kind=r4_kind) :: r4_value
    real(kind=r8_kind) :: r8_value
    integer(kind=i4_kind) :: i4_value
    integer(kind=i8_kind) :: i8_value

    call unpack_string(name)
    call unpack_words(value_type)
    select case (value_type(1))
    case (ATTRIBUTE_R4)
      r4_value = transfer(unpacked_words(r4_value), r4_value)
      call register_global_attribute(fileobj, trim(name), r4_value)
    case (ATTRIBUTE_R8)
      r8_value = transfer(unpacked_words(r8_value), r8_value)
      call register_global_attribute(fileobj, trim(name), r8_value)
    case (ATTRIBUTE_I4)
      i4_value = transfer(unpacked_words(i4_value), i4_value)
      call register_global_attribute(fileobj, trim(name), i4_value)
    case (ATTRIBUTE_I8)
      i8_value = transfer(unpacked_words(i8_value), i8_value)
      call register_global_attribute(fileobj, trim(name), i8_value)
    case (ATTRIBUTE_STR)
      call unpack_words(str_len)
      call unpack_string(str_value)
      call register_global_attribute(fileobj, trim(name), str_value, str_len=str_len(1))
    case default
      call mpp_error(FATAL, "offloading: unknown global attribute type received")
    end select
  end subroutine

  subroutine require_setup()
    if (.not. allocated(all_pes)) call mpp_error(FATAL, "offloading: call setup_offloading_pelists first")
  end subroutine

  subroutine require_model()
    call require_setup()
    if (.not. any(model_pes == mpp_pe())) call mpp_error(FATAL, "offloading: producer must run on model PEs")
    call mpp_set_current_pelist(model_pes, no_sync=.true.)
  end subroutine

  integer function require_file(fileobj) result(file_id)
    type(FmsNetcdfDomainFile_t), intent(in) :: fileobj

    call require_model()
    file_id = fileobj%offloading_obj_in%id
    call check_file_id(file_id)
    if (.not. file_is_open(file_id)) call mpp_error(FATAL, "offloading: no open file")
  end function

  integer function find_free_file_id() result(file_id)
    integer :: i

    file_id = 0
    do i = 1, max_files
      if (.not. file_is_open(i)) then
        file_id = i
        return
      endif
    enddo
  end function

  subroutine check_file_id(file_id)
    integer, intent(in) :: file_id

    if (file_id < 1 .or. file_id > max_files) call mpp_error(FATAL, "offloading: invalid file id")
  end subroutine

  subroutine require_writer_file(file_id)
    integer, intent(in) :: file_id

    call check_file_id(file_id)
    if (.not. file_is_open(file_id)) call mpp_error(FATAL, "offloading: operation received before file open")
  end subroutine

  subroutine check_string(value)
    character(len=*), intent(in) :: value
    if (len_trim(value) == 0 .or. len_trim(value) > string_length) &
      call mpp_error(FATAL, "offloading: names must contain 1 to 255 characters")
  end subroutine

  !> Called by all model PEs; only the model root builds and sends messages.
  subroutine start_message(operation)
    integer, intent(in) :: operation
    logical :: all_done
    integer :: k

    ! Let MPI progress field data sends that are still in flight.
    if (allocated(send_requests)) then
      do k = 1, size(send_requests, 2)
        call MPI_Testall(size(dest_ranks), send_requests(:,k), all_done, MPI_STATUSES_IGNORE)
      enddo
    endif
    message = 0
    message(1) = operation
    message_position = 1
  end subroutine

  !> Send the current message from the model root to every writer without waiting for delivery.
  !! Only waits if control_slots messages to a writer are still in flight.
  subroutine post_message()
    integer :: k, w

    if (mpp_pe() /= model_pes(1)) return
    k = next_control
    call MPI_Waitall(size(writer_pes), control_requests(:,k), MPI_STATUSES_IGNORE)
    next_control = mod(next_control, control_slots) + 1
    control_buffers(:,k) = message
    do w = 1, size(writer_pes)
      call MPI_Isend(control_buffers(:,k), message_position, MPI_INTEGER, writer_ranks(w), protocol_tag, &
                     offload_comm, control_requests(w,k))
    enddo
  end subroutine

  !> Writers: receive the next message from the model root and return its operation.
  integer function recv_message()
    call MPI_Recv(message, control_length, MPI_INTEGER, model_root_rank, protocol_tag, offload_comm, &
                  MPI_STATUS_IGNORE)
    recv_message = message(1)
    message_position = 1
  end function

  subroutine pack_words(words)
    integer, intent(in) :: words(:)
    message(message_position+1:message_position+size(words)) = words
    message_position = message_position + size(words)
  end subroutine

  subroutine pack_scalar(value)
    class(*), intent(in) :: value
    integer :: words(storage_size(value)/storage_size(0))

    words = transfer(value, words)
    call pack_words(words)
  end subroutine

  subroutine unpack_words(words)
    integer, intent(out) :: words(:)
    words = message(message_position+1:message_position+size(words))
    message_position = message_position + size(words)
  end subroutine

  !> Next words of the message, as many as needed to hold mold.
  function unpacked_words(mold) result(words)
    class(*), intent(in) :: mold
    integer, allocatable :: words(:)
    allocate(words(storage_size(mold)/storage_size(0)))
    call unpack_words(words)
  end function

  subroutine pack_string(value)
    character(len=*), intent(in) :: value
    character(len=string_words*storage_size(0)/8) :: padded

    padded = value
    call pack_words(transfer(padded, (/0/), string_words))
  end subroutine

  subroutine unpack_string(value)
    character(len=string_length), intent(out) :: value
    integer :: words(string_words)
    character(len=string_words*storage_size(0)/8) :: padded

    call unpack_words(words)
    padded = transfer(words, padded)
    value = padded
  end subroutine

  !> Rank of a PE in the offload communicator (pelists are in increasing order).
  integer function rank_of(pe)
    integer, intent(in) :: pe
    rank_of = findloc(all_pes, pe, 1) - 1
  end function
#endif
end module offloading_io_mod
