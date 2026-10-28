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
!! Numeric data may have up to four dimensions, plus an optional unlimited record dimension.
!! Supported kinds are r4, r8, i4 and i8; character data may be scalar or one-dimensional.
!! Scalar/vector data are replicated from the model root. Higher-rank data are decomposed in x/y.
!! As in fms2_io, int64 storage requires the netcdf4 file format.
module offloading_io_mod
#ifdef use_libMPI
!! TODO mpi calls should all go through mpp and have a non f08 interface ideally
  use mpi_f08, only: MPI_Comm, MPI_Request, MPI_Datatype, MPI_Isend, MPI_Irecv, MPI_Recv, MPI_Send, MPI_Wait, &
                     MPI_Waitall, MPI_Testall, MPI_Type_create_subarray, MPI_Type_commit, MPI_Type_free, &
                     MPI_INTEGER, MPI_INTEGER4, MPI_INTEGER8, MPI_CHARACTER, MPI_REAL4, MPI_REAL8, &
                     MPI_ORDER_FORTRAN, MPI_STATUS_IGNORE, MPI_DATATYPE_NULL, &
                     MPI_STATUSES_IGNORE, MPI_REQUEST_NULL
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_get_current_pelist, mpp_declare_pelist, mpp_set_current_pelist, &
                    mpp_sync, mpp_sum, mpp_error, FATAL, input_nml_file
  use fms_mod, only: check_nml_error
  use mpp_domains_mod
  use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, register_axis, register_field, write_data, &
                         register_global_attribute, unlimited
  use platform_mod, only: r4_kind, r8_kind, i4_kind, i8_kind

  implicit none
  private

  !! namelist variables
  !> Max number of io messages that can be in flight to each writer
  integer :: max_concurrent_messages = 64
  !> Maximum number of offloaded files that can be open simultaneously.
  integer :: max_files = 10
  namelist /offloading_io_nml/ max_files, max_concurrent_messages

  !> opcode enumerations to send the writers
  integer, parameter :: IO_OP_OPEN = 1, IO_OP_REGISTER_AXIS = 2, IO_OP_REGISTER_FIELD = 3
  integer, parameter :: IO_OP_WRITE = 4, IO_OP_CLOSE = 5, IO_OP_GLOBAL_ATTRIBUTE = 6, IO_OP_SHUTDOWN = 7
  integer, parameter :: ATTRIBUTE_R4 = 1, ATTRIBUTE_R8 = 2, ATTRIBUTE_I4 = 3, ATTRIBUTE_I8 = 4, ATTRIBUTE_STR = 5
  integer, parameter :: FIELD_R4 = 1, FIELD_R8 = 2, FIELD_I4 = 3, FIELD_I8 = 4, FIELD_STR = 5
  integer, parameter :: numeric_types = 4
  type(MPI_Datatype) :: field_mpi_types(numeric_types)
  !> These are arbitrary values used to tag the data flow in mpi communication
  !! protocol -> sending the opcodes, data -> sending the actual data to write
  integer, parameter :: protocol_tag = 24001, data_tag = 24002, string_length = 255
  integer, parameter :: default_max_pending_writes = 4
  !> Strings are packed into integers, 4 characters per word for 4 byte integers.
  integer, parameter :: string_words = (string_length + storage_size(0)/8) / (storage_size(0)/8)
  !> A control message holds a field name, storage type, up to five axis names and integer words.
  integer, parameter :: control_length = 7*string_words + 16
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
  type :: field_buffer_t
    class(*), allocatable :: values(:,:,:,:)
  end type
  type(field_buffer_t), allocatable, asynchronous :: send_buffers(:)
  type(MPI_Request), allocatable, asynchronous :: send_requests(:,:)

  ! Writer PE state: its part of the tile, the model ranks that overlap it and where their parts go.
  type(domain2d), save :: writer_domain
  type(FmsNetcdfDomainFile_t), allocatable :: writer_fileobjs(:)

  !> Private type to store information about the axes for each writer's files
  type :: writer_axes_t
    integer, allocatable :: directions(:)
    integer, allocatable :: sizes(:)
    character(len=string_length), allocatable :: names(:)
  end type
  type(writer_axes_t), allocatable :: writer_axes(:)
  integer, allocatable :: source_ranks(:)
  type(MPI_Datatype), allocatable :: recv_types(:,:)
  type(MPI_Request), allocatable, asynchronous :: recv_requests(:)
  integer :: writer_block_size(2) = 0
  class(*), allocatable, asynchronous :: block_data(:,:,:,:)

  public :: offloading_io_init, send_domain_metadata, recv_domain_metadata, run_offload_writer
  public :: open_offloaded_file, close_offloaded_file, offloading_io_end
  public :: register_axis_offloaded_file, register_axes_offloaded_file
  public :: register_field_offloaded_file, write_field_offloaded_file, global_metadata_offload

  !> Write r4, r8, i4 and i8 data of rank zero through four, or scalar/1-D character data.
  interface write_field_offloaded_file
    module procedure write_field_offloaded_file_0d
    module procedure write_field_offloaded_file_1d
    module procedure write_field_offloaded_file_2d
    module procedure write_field_offloaded_file_3d
    module procedure write_field_offloaded_file_4d
  end interface

contains

  !> Collective on the initial pelist; reserve its last nwriters PEs for the writers.
  !! nwriters must be a multiple of the number of tiles in the model domain (6 for a cubed sphere).
  subroutine offloading_io_init(global_pelist, model_pelist, writer_pelist, nwriters)
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
    field_mpi_types = [MPI_REAL4, MPI_REAL8, MPI_INTEGER4, MPI_INTEGER8]
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
    allocate(send_buffers(nslots))
    allocate(send_requests(size(writer_pes), nslots))
    send_requests = MPI_REQUEST_NULL
    next_send = 1
    if (mpp_pe() == model_pes(1)) then
      allocate(control_buffers(control_length, max_concurrent_messages))
      allocate(control_requests(size(writer_pes), max_concurrent_messages))
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
        allocate(writer_fileobjs(max_files), writer_axes(max_files))
    global_domain_size = global_size
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    metadata_ready = .true.
  end subroutine

  !> Collective on all PEs. Shares the model PEs' tiles and compute bounds, defines the writers'
  !! decomposition of the tiles, and sets up the MPI types that move each overlap between a model
  !! PE and a writer.
  subroutine exchange_layout(info)
    integer, intent(in) :: info(layout_info_size) !< this PE's layout information (zeros on writers)
    integer, allocatable :: table(:,:), writer_table(:,:) !< holds domain information for the model and writer PEs
                                                          !! first index is layout info element, second is PE 
                                                          !! 1 = tile id, 2 = is, 3 = ie, 4 = js, 5 = je,
                                                          !! 6 = global nx, 7 = global ny, 8 = ntiles,
                                                          !! 9-10 = writer layout, 11-12 = writer io layout
    integer, allocatable :: global_indices(:,:), layouts(:,:)
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
      allocate(recv_requests(max(1, size(source_ranks))))
      recv_requests = MPI_REQUEST_NULL
      writer_block_size = [writer_table(3,my_writer)-writer_table(2,my_writer)+1, &
               writer_table(5,my_writer)-writer_table(4,my_writer)+1]
      ncells = product(writer_block_size)
    endif
    if (overlap_cells /= ncells) call mpp_error(FATAL, "offloading: model and writer domains do not match")

  contains

    !> Find the PEs in others (same tile) whose bounds overlap, and build a subarray type for each
    !! overlap, relative to my own bounds.
    subroutine find_overlaps(mine, others, other_pes, ranks, types)
      integer, intent(in) :: mine(:), others(:,:), other_pes(:)
      integer, allocatable, intent(out) :: ranks(:)
      type(MPI_Datatype), allocatable, intent(out) :: types(:,:)
      integer :: m, k, field_kind

      allocate(ranks(0), types(numeric_types,0))
      overlap_cells = 0
      do m = 1, size(other_pes)
        if (others(1,m) /= mine(1)) cycle
        overlap = (/max(mine(2), others(2,m)), min(mine(3), others(3,m)), &
                    max(mine(4), others(4,m)), min(mine(5), others(5,m))/)
        if (overlap(1) > overlap(2) .or. overlap(3) > overlap(4)) cycle
        k = size(ranks) + 1
        ranks = (/ranks, rank_of(other_pes(m))/)
        types = reshape([types, (MPI_DATATYPE_NULL, field_kind=1,numeric_types)], [numeric_types, k])
        do field_kind = 1, numeric_types
          call MPI_Type_create_subarray(2, [mine(3)-mine(2)+1, mine(5)-mine(4)+1], &
                                        [overlap(2)-overlap(1)+1, overlap(4)-overlap(3)+1], &
                                        [overlap(1)-mine(2), overlap(3)-mine(4)], MPI_ORDER_FORTRAN, &
                                        field_mpi_types(field_kind), types(field_kind,k))
          call MPI_Type_commit(types(field_kind,k))
        enddo
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

  !> For x/y, axis_size is the tile size. Omit cart_name or pass a blank for other axes;
  !! axis_size is then a fixed length or fms2_io_mod's unlimited value.
  subroutine register_axis_offloaded_file(axis_name, axis_size, fileobj, cart_name)
    character(len=*), intent(in) :: axis_name
    integer, intent(in) :: axis_size
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    character(len=*), intent(in), optional :: cart_name
    integer :: direction, file_id
    character(len=1) :: cart

    file_id = require_file(fileobj)
    call check_string(axis_name)
    cart = ""
    if (present(cart_name)) then
      if (len_trim(cart_name) > 1) call mpp_error(FATAL, "offloading: axis cart_name must be x, y or blank")
      cart = cart_name
    endif
    direction = 0
    select case (cart)
    case ("x")
      direction = 1
    case ("y")
      direction = 2
    case ("")
    case default
      call mpp_error(FATAL, "offloading: axis cart_name must be x, y or blank")
    end select
    if (direction /= 0) then
      if (axis_size /= global_size(direction)) call mpp_error(FATAL, "offloading: axis size differs from domain")
    endif
    call start_message(IO_OP_REGISTER_AXIS)
    call pack_words((/file_id/))
    call pack_string(axis_name)
    call pack_string(cart)
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
      call register_axis_offloaded_file(axis_names(i), axis_sizes(i), fileobj, cart_names(i))
    enddo
  end subroutine

  !> axis_names identifies previously registered axes. Decomposed fields use x, y, extra axes,
  !! then optionally an unlimited axis; scalar/vector fields use their own registered axes.
  !! field_type defaults to double. Omitting axis_names registers a scalar field.
  subroutine register_field_offloaded_file(field_name, fileobj, field_type, axis_names)
    character(len=*), intent(in) :: field_name
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    character(len=*), intent(in), optional :: field_type
    character(len=*), intent(in), optional :: axis_names(:)
    integer :: file_id, ndims, dims
    character(len=string_length) :: storage_type

    file_id = require_file(fileobj)
    call check_string(field_name)
    ndims = 0
    if (present(axis_names)) then
      ndims = size(axis_names)
      if (ndims > 5) call mpp_error(FATAL, "offloading: fields support at most five dimensions")
      do dims = 1, ndims
        call check_string(axis_names(dims))
      enddo
    endif
    storage_type = "double"
    if (present(field_type)) storage_type = field_type
    call check_string(storage_type)
    call start_message(IO_OP_REGISTER_FIELD)
    call pack_words((/file_id/))
    call pack_string(field_name)
    call pack_words([ndims])
    do dims = 1, ndims
      call pack_string(axis_names(dims))
    enddo
    call pack_string(storage_type)
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

  !> Scalar and vector data are replicated from the model root to every writer.
  subroutine write_field_offloaded_file_0d(field_name, field_data, domain, fileobj, unlim_dim_level)
    character(len=*), intent(in) :: field_name
    class(*), intent(in) :: field_data
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer, intent(in), optional :: unlim_dim_level

    select type (field_data)
    type is (real(kind=r4_kind))
      call post_field(field_name, reshape([field_data], [1,1,1,1]), 0, domain, fileobj, unlim_dim_level)
    type is (real(kind=r8_kind))
      call post_field(field_name, reshape([field_data], [1,1,1,1]), 0, domain, fileobj, unlim_dim_level)
    type is (integer(kind=i4_kind))
      call post_field(field_name, reshape([field_data], [1,1,1,1]), 0, domain, fileobj, unlim_dim_level)
    type is (integer(kind=i8_kind))
      call post_field(field_name, reshape([field_data], [1,1,1,1]), 0, domain, fileobj, unlim_dim_level)
    type is (character(len=*))
      call post_field(field_name, reshape([field_data], [1,1,1,1]), 0, domain, fileobj, unlim_dim_level)
    class default
      call mpp_error(FATAL, "offloading: unsupported field data type")
    end select
  end subroutine

  subroutine write_field_offloaded_file_1d(field_name, field_data, domain, fileobj, unlim_dim_level)
    character(len=*), intent(in) :: field_name
    class(*), contiguous, target, intent(in) :: field_data(:)
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer, intent(in), optional :: unlim_dim_level
    class(*), pointer :: field4(:,:,:,:)

    field4(1:size(field_data),1:1,1:1,1:1) => field_data
    call post_field(field_name, field4, 1, domain, fileobj, unlim_dim_level)
  end subroutine

  !> Copy compute-domain data into a send ring. Only x/y are decomposed; halos are not sent.
  !! unlim_dim_level selects a positive record index; omission preserves write_data's default.
  subroutine write_field_offloaded_file_2d(field_name, field_data, domain, fileobj, unlim_dim_level)
    character(len=*), intent(in) :: field_name
    class(*), contiguous, target, intent(in) :: field_data(:,:)
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer, intent(in), optional :: unlim_dim_level
    class(*), pointer :: field4(:,:,:,:)

    field4(1:size(field_data,1),1:size(field_data,2),1:1,1:1) => field_data
    call post_field(field_name, field4, 2, domain, fileobj, unlim_dim_level)
  end subroutine

  subroutine write_field_offloaded_file_3d(field_name, field_data, domain, fileobj, unlim_dim_level)
    character(len=*), intent(in) :: field_name
    class(*), contiguous, target, intent(in) :: field_data(:,:,:)
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer, intent(in), optional :: unlim_dim_level
    class(*), pointer :: field4(:,:,:,:)

    field4(1:size(field_data,1),1:size(field_data,2),1:size(field_data,3),1:1) => field_data
    call post_field(field_name, field4, 3, domain, fileobj, unlim_dim_level)
  end subroutine

  subroutine write_field_offloaded_file_4d(field_name, field_data, domain, fileobj, unlim_dim_level)
    character(len=*), intent(in) :: field_name
    class(*), intent(in) :: field_data(:,:,:,:)
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer, intent(in), optional :: unlim_dim_level

    call post_field(field_name, field_data, 4, domain, fileobj, unlim_dim_level)
  end subroutine

  subroutine post_field(field_name, field_data, data_rank, domain, fileobj, unlim_dim_level)
    character(len=*), intent(in) :: field_name
    class(*), intent(in) :: field_data(:,:,:,:)
    integer, intent(in) :: data_rank
    type(domain2d), intent(inout) :: domain
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    integer, intent(in), optional :: unlim_dim_level
    integer :: file_id, level, field_kind, char_length, ioff, joff, nx, ny, slot, dest, count, destination
    logical :: replicated
    type(MPI_Datatype) :: datatype

    file_id = require_file(fileobj)
    call check_string(field_name)
    level = 0
    if (present(unlim_dim_level)) then
      if (unlim_dim_level < 1) call mpp_error(FATAL, "offloading: unlim_dim_level must be positive")
      level = unlim_dim_level
    endif
    char_length = 0
    select type (field_data)
    type is (real(kind=r4_kind))
      field_kind = FIELD_R4
    type is (real(kind=r8_kind))
      field_kind = FIELD_R8
    type is (integer(kind=i4_kind))
      field_kind = FIELD_I4
    type is (integer(kind=i8_kind))
      field_kind = FIELD_I8
    type is (character(len=*))
      if (data_rank > 1) call mpp_error(FATAL, "offloading: character data must be scalar or one-dimensional")
      if (level /= 0) call mpp_error(FATAL, "offloading: character writes do not support unlimited levels")
      field_kind = FIELD_STR
      char_length = len(field_data)
    class default
      call mpp_error(FATAL, "offloading: unsupported field data type")
    end select
    replicated = data_rank < 2
    if (replicated .and. mpp_pe() /= model_pes(1)) return
    ioff = 0
    joff = 0
    nx = size(field_data,1)
    ny = size(field_data,2)
    if (.not. replicated) call get_field_offsets([nx,ny], domain, ioff, joff, nx, ny)

    slot = next_send
    call MPI_Waitall(size(send_requests,1), send_requests(:,slot), MPI_STATUSES_IGNORE)
    next_send = mod(next_send, size(send_requests, 2)) + 1
    if (allocated(send_buffers(slot)%values)) deallocate(send_buffers(slot)%values)
    allocate(send_buffers(slot)%values, source=field_data(ioff+1:ioff+nx,joff+1:joff+ny,:,:))
    count = size(field_data,3)*size(field_data,4)
    if (replicated) count = size(field_data)
    if (field_kind == FIELD_STR) count = count*char_length
    do dest = 1, merge(size(writer_ranks), size(dest_ranks), replicated)
      if (replicated) then
        destination = writer_ranks(dest)
        datatype = MPI_CHARACTER
        if (field_kind /= FIELD_STR) datatype = field_mpi_types(field_kind)
      else
        destination = dest_ranks(dest)
        datatype = send_types(field_kind,dest)
      endif
      select type (buffer => send_buffers(slot)%values)
      type is (real(kind=r4_kind))
        call MPI_Isend(buffer, count, datatype, destination, data_tag, offload_comm, send_requests(dest,slot))
      type is (real(kind=r8_kind))
        call MPI_Isend(buffer, count, datatype, destination, data_tag, offload_comm, send_requests(dest,slot))
      type is (integer(kind=i4_kind))
        call MPI_Isend(buffer, count, datatype, destination, data_tag, offload_comm, send_requests(dest,slot))
      type is (integer(kind=i8_kind))
        call MPI_Isend(buffer, count, datatype, destination, data_tag, offload_comm, send_requests(dest,slot))
      type is (character(len=*))
        call MPI_Isend(buffer, count, datatype, destination, data_tag, offload_comm, send_requests(dest,slot))
      end select
    enddo
    call start_message(IO_OP_WRITE)
    call pack_words((/file_id/))
    call pack_string(field_name)
    call pack_words([field_kind, level, data_rank, nx, ny, size(field_data,3), size(field_data,4), char_length])
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
      call MPI_Waitall(size(send_requests,1), send_requests(:,w), MPI_STATUSES_IGNORE)
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
      do w = 1, max_concurrent_messages
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
      do d = 1, numeric_types
        call MPI_Type_free(send_types(d,w))
      enddo
    enddo
    deallocate(send_buffers, send_requests, dest_ranks, send_types)
    deallocate(file_is_open)
    session_finished = .true.
  end subroutine

  !> Enters the main write loop for offloaded writer PE's to receive data and control messages from the model PEs
  !! To be called by every writer after recv_domain_metadata.
  subroutine run_offload_writer()
    integer :: operation, file_id, file_id_words(1), axis_size(1), direction, field_dims(5), field_kind(1), n
    integer :: nfield_dims(1), level(1), data_info(6), field_kind_index
    character(len=string_length) :: name, cart, dimensions(5)

    ! checks for proper setup (offloading_io_init,send/recv_domain_metadata)
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
        allocate(writer_axes(file_id)%directions(0), writer_axes(file_id)%sizes(0), writer_axes(file_id)%names(0))
      case (IO_OP_REGISTER_AXIS)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call unpack_string(name)
        call unpack_string(cart)
        call unpack_words(axis_size)
        select case (trim(cart))
        case ("x")
          direction = 1
        case ("y")
          direction = 2
        case ("")
          direction = 0
        case default
          call mpp_error(FATAL, "offloading: invalid axis direction received")
        end select
        if (direction /= 0) then
          if (axis_size(1) /= global_size(direction)) &
            call mpp_error(FATAL, "offloading: invalid axis size received")
          if (any(writer_axes(file_id)%directions == direction)) &
            call mpp_error(FATAL, "offloading: duplicate axis direction")
          call register_axis(writer_fileobjs(file_id), trim(name), trim(cart))
        else
          call register_axis(writer_fileobjs(file_id), trim(name), axis_size(1))
        endif
        writer_axes(file_id)%names = [writer_axes(file_id)%names, name]
        writer_axes(file_id)%directions = [writer_axes(file_id)%directions, direction]
        writer_axes(file_id)%sizes = [writer_axes(file_id)%sizes, axis_size(1)]
      case (IO_OP_REGISTER_FIELD)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call unpack_string(name)
        call unpack_words(nfield_dims)
        if (nfield_dims(1) < 0 .or. nfield_dims(1) > 5) &
          call mpp_error(FATAL, "offloading: invalid field dimension count received")
        do n = 1, nfield_dims(1)
          call unpack_string(dimensions(n))
          field_dims(n) = findloc(writer_axes(file_id)%names, dimensions(n), dim=1)
          if (field_dims(n) == 0) &
            call mpp_error(FATAL, "offloading: field refers to unregistered axis "//trim(dimensions(n)))
        enddo
        call unpack_string(cart)
        if (nfield_dims(1) >= 2) then
          if (any(writer_axes(file_id)%directions(field_dims(:nfield_dims(1))) /= 0)) then
            if (any(writer_axes(file_id)%directions(field_dims(:2)) /= [1,2])) &
              call mpp_error(FATAL, "offloading: decomposed field dimensions must describe x then y")
            if (any(writer_axes(file_id)%directions(field_dims(3:nfield_dims(1))) /= 0)) &
              call mpp_error(FATAL, "offloading: only the first two field dimensions may be decomposed")
          endif
        endif
        if (any(writer_axes(file_id)%sizes(field_dims(:nfield_dims(1)-1)) == unlimited)) &
          call mpp_error(FATAL, "offloading: unlimited dimension must be last")
        if (nfield_dims(1) == 5) then
          if (writer_axes(file_id)%sizes(field_dims(5)) /= unlimited) &
            call mpp_error(FATAL, "offloading: fifth field dimension must be unlimited")
        endif
        if (nfield_dims(1) == 0) then
          call register_field(writer_fileobjs(file_id), trim(name), trim(cart))
        else
          call register_field(writer_fileobjs(file_id), trim(name), trim(cart), dimensions(:nfield_dims(1)))
        endif
      case (IO_OP_WRITE)
        call unpack_words(file_id_words)
        file_id = file_id_words(1)
        call require_writer_file(file_id)
        call unpack_string(name)
        call unpack_words(field_kind)
        call unpack_words(level)
        call unpack_words(data_info)
        call receive_field(field_kind(1), data_info)
        if (level(1) > 0) then
          call write_received_field(writer_fileobjs(file_id), trim(name), data_info(1), level(1))
        else
          call write_received_field(writer_fileobjs(file_id), trim(name), data_info(1))
        endif
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
        deallocate(writer_axes(file_id)%directions, writer_axes(file_id)%sizes, writer_axes(file_id)%names)
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
      do field_kind_index = 1, numeric_types
        call MPI_Type_free(recv_types(field_kind_index,n))
      enddo
    enddo
    if (allocated(block_data)) deallocate(block_data)
    deallocate(source_ranks, recv_types, recv_requests)
    deallocate(writer_fileobjs, writer_axes)
    deallocate(file_is_open)
  end subroutine

  subroutine receive_field(field_kind, data_info)
    integer, intent(in) :: field_kind, data_info(6)
    integer :: dims(4), count, source, source_rank, nsources
    logical :: replicated
    type(MPI_Datatype) :: datatype

    if (data_info(1) < 0 .or. data_info(1) > 4) call mpp_error(FATAL, "offloading: invalid data rank received")
    dims = data_info(2:5)
    replicated = data_info(1) < 2
    if (.not. replicated) dims(:2) = writer_block_size
    if (any(dims < 1)) call mpp_error(FATAL, "offloading: empty field dimensions are not supported")
    if (allocated(block_data)) deallocate(block_data)
    select case (field_kind)
    case (FIELD_R4)
      allocate(real(kind=r4_kind) :: block_data(dims(1),dims(2),dims(3),dims(4)))
    case (FIELD_R8)
      allocate(real(kind=r8_kind) :: block_data(dims(1),dims(2),dims(3),dims(4)))
    case (FIELD_I4)
      allocate(integer(kind=i4_kind) :: block_data(dims(1),dims(2),dims(3),dims(4)))
    case (FIELD_I8)
      allocate(integer(kind=i8_kind) :: block_data(dims(1),dims(2),dims(3),dims(4)))
    case (FIELD_STR)
      if (.not. replicated) call mpp_error(FATAL, "offloading: invalid character data rank received")
      allocate(character(len=data_info(6)) :: block_data(dims(1),dims(2),dims(3),dims(4)))
    case default
      call mpp_error(FATAL, "offloading: invalid field type received")
    end select
    count = dims(3)*dims(4)
    nsources = size(source_ranks)
    if (replicated) then
      count = product(dims)
      if (field_kind == FIELD_STR) count = count*data_info(6)
      nsources = 1
    endif
    do source = 1, nsources
      if (replicated) then
        source_rank = model_root_rank
        datatype = MPI_CHARACTER
        if (field_kind /= FIELD_STR) datatype = field_mpi_types(field_kind)
      else
        source_rank = source_ranks(source)
        datatype = recv_types(field_kind,source)
      endif
      select type (buffer => block_data)
      type is (real(kind=r4_kind))
        call MPI_Irecv(buffer, count, datatype, source_rank, data_tag, offload_comm, recv_requests(source))
      type is (real(kind=r8_kind))
        call MPI_Irecv(buffer, count, datatype, source_rank, data_tag, offload_comm, recv_requests(source))
      type is (integer(kind=i4_kind))
        call MPI_Irecv(buffer, count, datatype, source_rank, data_tag, offload_comm, recv_requests(source))
      type is (integer(kind=i8_kind))
        call MPI_Irecv(buffer, count, datatype, source_rank, data_tag, offload_comm, recv_requests(source))
      type is (character(len=*))
        call MPI_Irecv(buffer, count, datatype, source_rank, data_tag, offload_comm, recv_requests(source))
      end select
    enddo
    call MPI_Waitall(size(recv_requests), recv_requests, MPI_STATUSES_IGNORE)
  end subroutine

  subroutine write_received_field(fileobj, field_name, data_rank, unlim_dim_level)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: data_rank
    integer, intent(in), optional :: unlim_dim_level

    select case (data_rank)
    case (0)
      call write_data(fileobj, field_name, block_data(1,1,1,1), unlim_dim_level=unlim_dim_level)
    case (1)
      call write_data(fileobj, field_name, block_data(:,1,1,1), unlim_dim_level=unlim_dim_level)
    case (2)
      call write_data(fileobj, field_name, block_data(:,:,1,1), unlim_dim_level=unlim_dim_level)
    case (3)
      call write_data(fileobj, field_name, block_data(:,:,:,1), unlim_dim_level=unlim_dim_level)
    case (4)
      call write_data(fileobj, field_name, block_data, unlim_dim_level=unlim_dim_level)
    end select
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
    if (.not. allocated(all_pes)) call mpp_error(FATAL, "offloading: call offloading_io_init first")
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
        call MPI_Testall(size(send_requests,1), send_requests(:,k), all_done, MPI_STATUSES_IGNORE)
      enddo
    endif
    message = 0
    message(1) = operation
    message_position = 1
  end subroutine

  !> Send the current message from the model root to every writer without waiting for delivery.
  !! Only waits if max_concurrent_messages messages to a writer are still in flight.
  subroutine post_message()
    integer :: k, w

    if (mpp_pe() /= model_pes(1)) return
    k = next_control
    call MPI_Waitall(size(writer_pes), control_requests(:,k), MPI_STATUSES_IGNORE)
    next_control = mod(next_control, max_concurrent_messages) + 1
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
