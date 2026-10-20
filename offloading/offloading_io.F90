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
!> @brief This module provides a method of offloading I/O operations to dedicated writer PEs.
!!
!! One writer PE is used per domain tile, so single-tile domains use one writer and cubed-sphere
!! domains (see mpp_define_cubic_domain) use six. Writer N writes tile N, and fms2_io adds the
!! ".tileN" suffix to the file name for multi-tile domains.
!!
!! The model side never waits on the writers during normal operation:
!! - Field data does not go through a model-side gather. Each model PE copies its compute domain
!!   into a buffer from a small ring and posts a nonblocking send directly to its tile's writer.
!!   A send only has to finish when its buffer is reused or the file is closed.
!! - Each operation (open, register, write, ...) is one packed control message from the model root
!!   to every writer, also sent without waiting from a ring of buffers.
!! Writers receive field data straight into the tile array through MPI subarray types, and each
!! tile is written by its own writer at the same time.
module offloading_io_mod
#ifdef use_libMPI
!! TODO mpi calls should all go through mpp and have a non f08 interface ideally
  use mpi_f08, only: MPI_Comm, MPI_Request, MPI_Datatype, MPI_Isend, MPI_Irecv, MPI_Recv, MPI_Send, MPI_Wait, &
                     MPI_Waitall, MPI_Testall, MPI_Type_create_subarray, MPI_Type_commit, MPI_Type_free, &
                     MPI_INTEGER, MPI_REAL4, MPI_REAL8, MPI_ORDER_FORTRAN, MPI_STATUS_IGNORE, &
                     MPI_STATUSES_IGNORE, MPI_REQUEST_NULL
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_get_current_pelist, mpp_declare_pelist, mpp_set_current_pelist, &
                    mpp_sync, mpp_sum, mpp_error, FATAL
  use mpp_domains_mod
  use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, register_axis, register_field, write_data, &
                         register_global_attribute
  use platform_mod, only: r4_kind, r8_kind, i4_kind, i8_kind

  implicit none
  private

  integer, parameter :: IO_OP_OPEN = 1, IO_OP_REGISTER_AXIS = 2, IO_OP_REGISTER_FIELD = 3
  integer, parameter :: IO_OP_WRITE = 4, IO_OP_CLOSE = 5, IO_OP_GLOBAL_ATTRIBUTE = 6
  integer, parameter :: ATTRIBUTE_R4 = 1, ATTRIBUTE_R8 = 2, ATTRIBUTE_I4 = 3, ATTRIBUTE_I8 = 4, ATTRIBUTE_STR = 5
  integer, parameter :: protocol_tag = 24001, data_tag = 24002, string_length = 255
  integer, parameter :: default_max_pending_writes = 4
  !> Strings are packed into integers, 4 characters per word for 4 byte integers.
  integer, parameter :: string_words = (string_length + storage_size(0)/8) / (storage_size(0)/8)
  !> A control message holds an operation and at most two strings plus a few integer words.
  integer, parameter :: control_length = 2*string_words + 8
  !> Number of control messages the model root can have in flight to each writer.
  integer, parameter :: control_slots = 64
  !> Columns of the layout table: tile, isc, iec, jsc, jec, tile x size, tile y size
  integer, parameter :: layout_info_size = 7
  integer, allocatable :: all_pes(:), model_pes(:), writer_pes(:)
  integer, allocatable :: writer_ranks(:)
  integer :: model_root_rank = -1
  integer :: global_size(2) = 0 !< size of each tile
  integer :: ntiles = 0
  logical :: metadata_ready = .false., file_is_open = .false., session_finished = .false.
  type(MPI_Comm) :: offload_comm
  type(MPI_Datatype) :: real_type

  ! Control message being packed (model root) or unpacked (writers).
  integer :: message(control_length), message_position = 0

  ! Model root state: ring of control messages in flight to the writers.
  integer :: next_control = 1
  integer, allocatable, asynchronous :: control_buffers(:,:)
  type(MPI_Request), allocatable, asynchronous :: control_requests(:,:)

  ! Model PE state: compute bounds, destination writer, and the ring of field data in flight.
  integer :: compute_bounds(4) = 0, writer_rank = -1, next_send = 1
  real, allocatable, asynchronous :: send_buffers(:,:,:)
  type(MPI_Request), allocatable, asynchronous :: send_requests(:)

  ! Writer PE state: the model ranks on this tile and where their blocks go in the tile array.
  integer, allocatable :: source_ranks(:)
  type(MPI_Datatype), allocatable :: recv_types(:)
  type(MPI_Request), allocatable, asynchronous :: recv_requests(:)
  real, allocatable, asynchronous :: tile_data(:,:)

  public :: setup_offloading_pelists, send_domain_metadata, recv_domain_metadata, run_offload_writer
  public :: open_offloaded_file, close_offloaded_file, close_offloaded_files
  public :: register_axis_offloaded_file, register_axes_offloaded_file
  public :: register_field_offloaded_file, write_field_offloaded_file, global_metadata_offload

contains

  !> Collective on the initial pelist; reserve its last nwriters PEs for the writers.
  !! nwriters must equal the number of tiles in the model domain (6 for a cubed sphere).
  subroutine setup_offloading_pelists(global_pelist, model_pelist, writer_pelist, nwriters)
    integer, allocatable, intent(out) :: global_pelist(:), model_pelist(:), writer_pelist(:)
    integer, intent(in), optional :: nwriters !< number of writer PEs, one per tile (default 1)
    integer :: n, comm_id, w

    if (allocated(all_pes)) call mpp_error(FATAL, "offloading: pelists already initialized")
    n = 1
    if (present(nwriters)) n = nwriters
    if (n < 1) call mpp_error(FATAL, "offloading: at least one writer PE is required")
    if (mpp_npes() <= n) call mpp_error(FATAL, "offloading: at least one model PE is required besides the writers")
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
    if (storage_size(1.0) == 64) then
      real_type = MPI_REAL8
    else
      real_type = MPI_REAL4
    endif
    global_pelist = all_pes
    model_pelist = model_pes
    writer_pelist = writer_pes
  end subroutine

  !> All model PEs call this once, paired with recv_domain_metadata on every writer.
  !! The domain must have one tile per writer PE, and each model PE must own exactly one tile.
  !! max_pending_writes sets how many field writes each model PE can have in flight (default 4).
  subroutine send_domain_metadata(model_domain, max_pending_writes)
    type(domain2d), intent(inout) :: model_domain
    integer, intent(in), optional :: max_pending_writes
    integer :: info(layout_info_size), tile_id(1), nslots

    call require_model()
    if (metadata_ready) call mpp_error(FATAL, "offloading: metadata already sent")
    if (mpp_get_ntile_count(model_domain) /= size(writer_pes)) &
      call mpp_error(FATAL, "offloading: the number of writer PEs must equal the number of domain tiles")
    if (mpp_get_current_ntile(model_domain) /= 1) &
      call mpp_error(FATAL, "offloading: each model PE must own exactly one tile")
    nslots = default_max_pending_writes
    if (present(max_pending_writes)) nslots = max_pending_writes
    if (nslots < 1) call mpp_error(FATAL, "offloading: max_pending_writes must be positive")
    tile_id = mpp_get_tile_id(model_domain)
    info(1) = tile_id(1)
    call mpp_get_compute_domain(model_domain, info(2), info(3), info(4), info(5))
    call mpp_get_global_domain(model_domain, xsize=info(6), ysize=info(7))
    call mpp_set_current_pelist(all_pes, no_sync=.true.)
    call mpp_broadcast_domain(model_domain)
    call exchange_layout(info)
    compute_bounds = info(2:5)
    writer_rank = writer_ranks(info(1))
    allocate(send_buffers(info(3)-info(2)+1, info(5)-info(4)+1, nslots))
    allocate(send_requests(nslots))
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
    global_domain_size = global_size
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    metadata_ready = .true.
  end subroutine

  !> Collective on all PEs. Shares each model PE's tile and compute bounds, checks that every tile
  !! is covered and has the same size, and sets up the receive types on the writers.
  subroutine exchange_layout(info)
    integer, intent(in) :: info(layout_info_size) !< this PE's layout information (zeros on writers)
    integer, allocatable :: table(:,:)
    integer :: n, k, my_tile, ncells

    ntiles = size(writer_pes)
    allocate(table(layout_info_size, size(model_pes)))
    table = 0
    do n = 1, size(model_pes)
      if (model_pes(n) == mpp_pe()) table(:,n) = info
    enddo
    call mpp_sum(table, size(table))
    global_size = table(6:7, 1)
    do n = 1, size(model_pes)
      if (table(1,n) < 1 .or. table(1,n) > ntiles) call mpp_error(FATAL, "offloading: invalid model tile id")
      if (any(table(6:7,n) /= global_size)) call mpp_error(FATAL, "offloading: all tiles must have the same size")
    enddo
    do n = 1, ntiles
      if (.not. any(table(1,:) == n)) call mpp_error(FATAL, "offloading: a tile has no model PEs")
    enddo

    my_tile = findloc(writer_pes, mpp_pe(), 1)
    if (my_tile == 0) return
    allocate(source_ranks(count(table(1,:) == my_tile)))
    allocate(recv_types(size(source_ranks)), recv_requests(size(source_ranks)))
    allocate(tile_data(global_size(1), global_size(2)))
    k = 0
    ncells = 0
    do n = 1, size(model_pes)
      if (table(1,n) /= my_tile) cycle
      k = k + 1
      source_ranks(k) = rank_of(model_pes(n))
      call MPI_Type_create_subarray(2, global_size, (/table(3,n)-table(2,n)+1, table(5,n)-table(4,n)+1/), &
                                    (/table(2,n)-1, table(4,n)-1/), MPI_ORDER_FORTRAN, real_type, recv_types(k))
      call MPI_Type_commit(recv_types(k))
      ncells = ncells + (table(3,n)-table(2,n)+1) * (table(5,n)-table(4,n)+1)
    enddo
    if (ncells /= product(global_size)) call mpp_error(FATAL, "offloading: compute domains do not cover the tile")
    recv_requests = MPI_REQUEST_NULL
  end subroutine

  subroutine open_offloaded_file(filename)
    character(len=*), intent(in) :: filename

    call require_model()
    if (.not. metadata_ready) call mpp_error(FATAL, "offloading: send domain metadata before opening a file")
    if (file_is_open .or. session_finished) call mpp_error(FATAL, "offloading: only one file per session is supported")
    call check_string(filename)
    call start_message(IO_OP_OPEN)
    call pack_string(filename)
    call post_message()
    file_is_open = .true.
  end subroutine

  !> axis_size is the size of one tile in the cart_name direction.
  subroutine register_axis_offloaded_file(axis_name, cart_name, axis_size)
    character(len=*), intent(in) :: axis_name, cart_name
    integer, intent(in) :: axis_size
    integer :: direction

    call require_file()
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
    call pack_string(axis_name)
    call pack_string(cart_name)
    call pack_words((/axis_size/))
    call post_message()
  end subroutine

  subroutine register_axes_offloaded_file(axis_names, cart_names, axis_sizes, domain)
    character(len=*), intent(in) :: axis_names(:), cart_names(:)
    integer, intent(in) :: axis_sizes(:)
    type(domain2d), intent(in) :: domain
    integer :: i, nx, ny

    call require_file()
    if (size(axis_names) /= size(cart_names) .or. size(axis_names) /= size(axis_sizes)) &
      call mpp_error(FATAL, "offloading: axis arrays must have equal sizes")
    call mpp_get_global_domain(domain, xsize=nx, ysize=ny)
    if (any(global_size /= (/nx, ny/))) call mpp_error(FATAL, "offloading: axis domain differs from metadata")
    do i = 1, size(axis_names)
      call register_axis_offloaded_file(axis_names(i), cart_names(i), axis_sizes(i))
    enddo
  end subroutine

  !> Dimension indices refer to registration order, and must describe x then y.
  subroutine register_field_offloaded_file(field_name, field_dims)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: field_dims(:)

    call require_file()
    call check_string(field_name)
    if (size(field_dims) /= 2) call mpp_error(FATAL, "offloading: only 2-D fields are supported")
    if (any(field_dims < 1) .or. any(field_dims > 2)) &
      call mpp_error(FATAL, "offloading: field dimension indices must be 1 or 2")
    call start_message(IO_OP_REGISTER_FIELD)
    call pack_string(field_name)
    call pack_words(field_dims)
    call post_message()
  end subroutine

  !> Register a scalar global attribute (r4, r8, i4, i8 or string) on the offloaded file.
  !! Only the model root's value is sent, and every tile file gets it.
  !! Can be called at any time while the file is open.
  subroutine global_metadata_offload(attribute_name, attribute_value)
    character(len=*), intent(in) :: attribute_name
    class(*), intent(in) :: attribute_value

    call require_file()
    call check_string(attribute_name)
    call start_message(IO_OP_GLOBAL_ATTRIBUTE)
    call pack_string(attribute_name)
    select type (attribute_value)
    type is (real(kind=r4_kind))
      call pack_words((/ATTRIBUTE_R4/))
      call pack_words(transfer(attribute_value, (/0/)))
    type is (real(kind=r8_kind))
      call pack_words((/ATTRIBUTE_R8/))
      call pack_words(transfer(attribute_value, (/0/)))
    type is (integer(kind=i4_kind))
      call pack_words((/ATTRIBUTE_I4/))
      call pack_words(transfer(attribute_value, (/0/)))
    type is (integer(kind=i8_kind))
      call pack_words((/ATTRIBUTE_I8/))
      call pack_words(transfer(attribute_value, (/0/)))
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
  subroutine write_field_offloaded_file(field_name, field_data, domain)
    character(len=*), intent(in) :: field_name
    real, intent(in) :: field_data(:,:)
    type(domain2d), intent(inout) :: domain
    integer :: nx, ny, isc, iec, jsc, jec, isd, ied, jsd, jed, ioff, joff, k

    call require_file()
    call check_string(field_name)
    call mpp_get_global_domain(domain, xsize=nx, ysize=ny)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)
    if (any(global_size /= (/nx, ny/)) .or. any(compute_bounds /= (/isc, iec, jsc, jec/))) &
      call mpp_error(FATAL, "offloading: field domain differs from metadata")
    if (all(shape(field_data) == (/iec-isc+1, jec-jsc+1/))) then
      ioff = 0
      joff = 0
    else if (all(shape(field_data) == (/ied-isd+1, jed-jsd+1/))) then
      ioff = isc - isd
      joff = jsc - jsd
    else
      call mpp_error(FATAL, "offloading: field data must cover the compute or data domain")
    endif

    ! Reuse the oldest buffer; this only waits if all max_pending_writes sends are still in flight.
    k = next_send
    call MPI_Wait(send_requests(k), MPI_STATUS_IGNORE)
    next_send = mod(next_send, size(send_requests)) + 1
    send_buffers(:,:,k) = field_data(ioff+1:ioff+iec-isc+1, joff+1:joff+jec-jsc+1)
    call MPI_Isend(send_buffers(:,:,k), size(send_buffers(:,:,k)), real_type, writer_rank, data_tag, &
                   offload_comm, send_requests(k))
    call start_message(IO_OP_WRITE)
    call pack_string(field_name)
    call post_message()
  end subroutine

  !> Close also stops the writers. Returns once every tile file is closed.
  subroutine close_offloaded_file()
    integer :: acknowledgement, w

    call require_file()
    call MPI_Waitall(size(send_requests), send_requests, MPI_STATUSES_IGNORE)
    call start_message(IO_OP_CLOSE)
    call post_message()
    if (mpp_pe() == model_pes(1)) then
      do w = 1, control_slots
        call MPI_Waitall(size(writer_pes), control_requests(:,w), MPI_STATUSES_IGNORE)
      enddo
      do w = 1, size(writer_pes)
        call MPI_Recv(acknowledgement, 1, MPI_INTEGER, writer_ranks(w), protocol_tag, offload_comm, &
                      MPI_STATUS_IGNORE)
        if (acknowledgement /= IO_OP_CLOSE) call mpp_error(FATAL, "offloading: invalid close acknowledgement")
      enddo
      deallocate(control_buffers, control_requests)
    endif
    call mpp_sync()
    deallocate(send_buffers, send_requests)
    file_is_open = .false.
    session_finished = .true.
  end subroutine

  subroutine close_offloaded_files()
    call close_offloaded_file()
  end subroutine

  !> Enters main write loop for an offloaded writer PE.
  !! To be called by every writer after recv_domain_metadata; writer N writes tile N.
  subroutine run_offload_writer()
    type(FmsNetcdfDomainFile_t) :: fileobj
    type(domain2d) :: domain
    integer :: operation, axis_size(1), naxes, direction, field_dims(2), axis_directions(2), n
    integer :: global_indices(4, size(writer_pes)), layout(2, size(writer_pes)), no_contact(0)
    character(len=string_length) :: name, cart, axis_names(2)

    call require_setup()
    if (.not. any(writer_pes == mpp_pe())) call mpp_error(FATAL, "offloading: consumer must run on a writer PE")
    if (.not. metadata_ready) call mpp_error(FATAL, "offloading: receive metadata before running writer")
    if (session_finished) call mpp_error(FATAL, "offloading: writer session already finished")
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    do n = 1, ntiles
      global_indices(:,n) = (/1, global_size(1), 1, global_size(2)/)
      layout(:,n) = (/1, 1/)
    enddo
    call mpp_define_mosaic(global_indices, layout, domain, ntiles, 0, no_contact, no_contact, no_contact, &
                           no_contact, no_contact, no_contact, no_contact, no_contact, no_contact, no_contact, &
                           writer_pes, writer_pes, name="offload writer mosaic")
    call mpp_define_io_domain(domain, (/1,1/))
    naxes = 0
    axis_directions = 0
    do
      operation = recv_message()
      if (operation /= IO_OP_OPEN .and. .not. file_is_open) &
        call mpp_error(FATAL, "offloading: operation received before file open")
      select case (operation)
      case (IO_OP_OPEN)
        if (file_is_open) call mpp_error(FATAL, "offloading: file already open")
        call unpack_string(name)
        if (.not. open_file(fileobj, trim(name), "overwrite", domain=domain)) &
          call mpp_error(FATAL, "offloading: unable to open "//trim(name))
        file_is_open = .true.
      case (IO_OP_REGISTER_AXIS)
        call unpack_string(name)
        call unpack_string(cart)
        call unpack_words(axis_size)
        if (naxes == 2) call mpp_error(FATAL, "offloading: only two axes are supported")
        select case (trim(cart))
        case ("x")
          direction = 1
        case ("y")
          direction = 2
        case default
          call mpp_error(FATAL, "offloading: invalid axis direction received")
        end select
        if (axis_size(1) /= global_size(direction)) call mpp_error(FATAL, "offloading: invalid axis size received")
        if (any(axis_directions(:naxes) == direction)) call mpp_error(FATAL, "offloading: duplicate axis direction")
        naxes = naxes + 1
        axis_names(naxes) = name
        axis_directions(naxes) = direction
        call register_axis(fileobj, trim(name), trim(cart))
      case (IO_OP_REGISTER_FIELD)
        call unpack_string(name)
        call unpack_words(field_dims)
        if (any(field_dims < 1) .or. any(field_dims > naxes)) &
          call mpp_error(FATAL, "offloading: field refers to unregistered axes")
        if (any(axis_directions(field_dims) /= (/1,2/))) &
          call mpp_error(FATAL, "offloading: field dimensions must describe x then y")
        call register_field(fileobj, trim(name), "double", axis_names(field_dims))
      case (IO_OP_WRITE)
        call unpack_string(name)
        do n = 1, size(source_ranks)
          call MPI_Irecv(tile_data, 1, recv_types(n), source_ranks(n), data_tag, offload_comm, recv_requests(n))
        enddo
        call MPI_Waitall(size(recv_requests), recv_requests, MPI_STATUSES_IGNORE)
        call write_data(fileobj, trim(name), tile_data)
      case (IO_OP_GLOBAL_ATTRIBUTE)
        call recv_global_attribute(fileobj)
      case (IO_OP_CLOSE)
        call close_file(fileobj)
        file_is_open = .false.
        session_finished = .true.
        call MPI_Send(IO_OP_CLOSE, 1, MPI_INTEGER, model_root_rank, protocol_tag, offload_comm)
        exit
      case default
        call mpp_error(FATAL, "offloading: unknown operation received")
      end select
    enddo
    call mpp_deallocate_domain(domain)
    do n = 1, size(recv_types)
      call MPI_Type_free(recv_types(n))
    enddo
    deallocate(source_ranks, recv_types, recv_requests, tile_data)
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

  subroutine require_file()
    call require_model()
    if (.not. file_is_open) call mpp_error(FATAL, "offloading: no open file")
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

    ! Let MPI progress field data sends that are still in flight.
    if (allocated(send_requests)) &
      call MPI_Testall(size(send_requests), send_requests, all_done, MPI_STATUSES_IGNORE)
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

  subroutine unpack_words(words)
    integer, intent(out) :: words(:)
    words = message(message_position+1:message_position+size(words))
    message_position = message_position + size(words)
  end subroutine

  !> Next words of the message, as many as needed to hold mold.
  function unpacked_words(mold) result(words)
    class(*), intent(in) :: mold
    integer :: words(storage_size(mold)/storage_size(0))
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
