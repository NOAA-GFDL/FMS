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
module offloading_io_mod
#ifdef use_libMPI
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_get_current_pelist, mpp_declare_pelist, mpp_set_current_pelist, &
                    mpp_broadcast, mpp_send, mpp_recv, mpp_sync_self, mpp_sync, mpp_error, FATAL
  use mpp_parameter_mod, only: GLOBAL_ROOT_ONLY, XUPDATE, YUPDATE
  use mpp_domains_mod
  use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, register_axis, register_field, write_data, &
                         register_global_attribute
  use platform_mod, only: r4_kind, r8_kind, i4_kind, i8_kind

  implicit none
  private

  integer, parameter :: IO_OP_OPEN = 1, IO_OP_REGISTER_AXIS = 2, IO_OP_REGISTER_FIELD = 3
  integer, parameter :: IO_OP_WRITE = 4, IO_OP_CLOSE = 5, IO_OP_GLOBAL_ATTRIBUTE = 6
  integer, parameter :: ATTRIBUTE_R4 = 1, ATTRIBUTE_R8 = 2, ATTRIBUTE_I4 = 3, ATTRIBUTE_I8 = 4, ATTRIBUTE_STR = 5
  integer, parameter :: protocol_tag = 24001, string_length = 255
  integer, allocatable :: all_pes(:), model_pes(:), writer_pes(:)
  integer :: global_size(2) = 0
  logical :: metadata_ready = .false., file_is_open = .false., session_finished = .false.

  public :: setup_offloading_pelists, send_domain_metadata, recv_domain_metadata, run_offload_writer
  public :: open_offloaded_file, close_offloaded_file, close_offloaded_files
  public :: register_axis_offloaded_file, register_axes_offloaded_file
  public :: register_field_offloaded_file, write_field_offloaded_file, global_metadata_offload

contains

  !> Collective on the initial pelist; reserve its last PE for the consumer.
  subroutine setup_offloading_pelists(global_pelist, model_pelist, writer_pelist)
    integer, allocatable, intent(out) :: global_pelist(:), model_pelist(:), writer_pelist(:)

    if (allocated(all_pes)) call mpp_error(FATAL, "offloading: pelists already initialized")
    if (mpp_npes() < 2) call mpp_error(FATAL, "offloading: at least one model PE and one writer PE are required")
    allocate(all_pes(mpp_npes()))
    call mpp_get_current_pelist(all_pes)
    model_pes = all_pes(:size(all_pes)-1)
    writer_pes = all_pes(size(all_pes):)
    call mpp_declare_pelist(all_pes)
    call mpp_declare_pelist(model_pes)
    call mpp_declare_pelist(writer_pes)
    global_pelist = all_pes
    model_pelist = model_pes
    writer_pelist = writer_pes
  end subroutine

  !> All model PEs call this once, paired with recv_domain_metadata on the writer.
  subroutine send_domain_metadata(model_domain)
    type(domain2d), intent(inout) :: model_domain

    call require_model()
    if (metadata_ready) call mpp_error(FATAL, "offloading: metadata already sent")
    if (mpp_get_ntile_count(model_domain) /= 1) &
      call mpp_error(FATAL, "offloading: only single-tile domains are supported")
    call mpp_get_global_domain(model_domain, xsize=global_size(1), ysize=global_size(2))
    call mpp_set_current_pelist(all_pes, no_sync=.true.)
    call mpp_broadcast_domain(model_domain)
    call mpp_broadcast(global_size, 2, model_pes(1))
    call mpp_set_current_pelist(model_pes, no_sync=.true.)
    metadata_ready = .true.
  end subroutine

  subroutine recv_domain_metadata(model_domain, global_domain_size)
    type(domain2d), intent(inout) :: model_domain
    integer, intent(out) :: global_domain_size(2)

    call require_setup()
    if (mpp_pe() /= writer_pes(1)) call mpp_error(FATAL, "offloading: receiver must be the writer PE")
    if (metadata_ready) call mpp_error(FATAL, "offloading: metadata already received")
    call mpp_set_current_pelist(all_pes, no_sync=.true.)
    call mpp_define_null_domain(model_domain)
    call mpp_broadcast_domain(model_domain)
    call mpp_broadcast(global_size, 2, model_pes(1))
    global_domain_size = global_size
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    metadata_ready = .true.
  end subroutine

  subroutine open_offloaded_file(filename)
    character(len=*), intent(in) :: filename

    call require_model()
    if (.not. metadata_ready) call mpp_error(FATAL, "offloading: send domain metadata before opening a file")
    if (file_is_open .or. session_finished) call mpp_error(FATAL, "offloading: only one file per session is supported")
    call check_string(filename)
    call begin_operation(IO_OP_OPEN)
    if (mpp_pe() == model_pes(1)) call send_string(filename)
    call end_operation()
    file_is_open = .true.
  end subroutine

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
    call begin_operation(IO_OP_REGISTER_AXIS)
    if (mpp_pe() == model_pes(1)) then
      call send_string(axis_name)
      call send_string(cart_name)
      call mpp_send(axis_size, writer_pes(1), tag=protocol_tag)
    endif
    call end_operation()
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
    call begin_operation(IO_OP_REGISTER_FIELD)
    if (mpp_pe() == model_pes(1)) then
      call send_string(field_name)
      call mpp_send(field_dims, 2, writer_pes(1), tag=protocol_tag)
    endif
    call end_operation()
  end subroutine

  !> Register a scalar global attribute (r4, r8, i4, i8 or string) on the offloaded file.
  !! Only the model root's value is sent. Can be called at any time while the file is open.
  subroutine global_metadata_offload(attribute_name, attribute_value)
    character(len=*), intent(in) :: attribute_name
    class(*), intent(in) :: attribute_value
    integer :: value_type, str_len

    call require_file()
    call check_string(attribute_name)
    str_len = 0
    select type (attribute_value)
    type is (real(kind=r4_kind))
      value_type = ATTRIBUTE_R4
    type is (real(kind=r8_kind))
      value_type = ATTRIBUTE_R8
    type is (integer(kind=i4_kind))
      value_type = ATTRIBUTE_I4
    type is (integer(kind=i8_kind))
      value_type = ATTRIBUTE_I8
    type is (character(len=*))
      value_type = ATTRIBUTE_STR
      str_len = len_trim(attribute_value)
      if (str_len > string_length) &
        call mpp_error(FATAL, "offloading: string attribute values must contain at most 255 characters")
    class default
      call mpp_error(FATAL, "offloading: unsupported global attribute type for "//trim(attribute_name))
    end select
    call begin_operation(IO_OP_GLOBAL_ATTRIBUTE)
    if (mpp_pe() == model_pes(1)) then
      call send_string(attribute_name)
      call mpp_send(value_type, writer_pes(1), tag=protocol_tag)
      select type (attribute_value)
      type is (real(kind=r4_kind))
        call mpp_send(attribute_value, writer_pes(1), tag=protocol_tag)
      type is (real(kind=r8_kind))
        call mpp_send(attribute_value, writer_pes(1), tag=protocol_tag)
      type is (integer(kind=i4_kind))
        call mpp_send(attribute_value, writer_pes(1), tag=protocol_tag)
      type is (integer(kind=i8_kind))
        call mpp_send(attribute_value, writer_pes(1), tag=protocol_tag)
      type is (character(len=*))
        call mpp_send(str_len, writer_pes(1), tag=protocol_tag)
        call send_string(attribute_value)
      end select
    endif
    call end_operation()
  end subroutine

  !> Collective on model PEs. FMS gathers compute cells, excluding any halos.
  subroutine write_field_offloaded_file(field_name, field_data, domain)
    character(len=*), intent(in) :: field_name
    real, intent(in) :: field_data(:,:)
    type(domain2d), intent(inout) :: domain
    real, allocatable :: global_data(:,:)
    integer :: nx, ny

    call require_file()
    call check_string(field_name)
    call mpp_get_global_domain(domain, xsize=nx, ysize=ny)
    if (any(global_size /= (/nx, ny/))) call mpp_error(FATAL, "offloading: field domain differs from metadata")
    allocate(global_data(nx, ny))
    call mpp_global_field(domain, field_data, global_data, flags=GLOBAL_ROOT_ONLY+XUPDATE+YUPDATE)
    call begin_operation(IO_OP_WRITE)
    if (mpp_pe() == model_pes(1)) then
      call send_string(field_name)
      call mpp_send(global_data, size(global_data), writer_pes(1), tag=protocol_tag)
    endif
    call end_operation()
  end subroutine

  !> Close also stops the consumer. The acknowledgement means the file is closed.
  subroutine close_offloaded_file()
    integer :: acknowledgement

    call require_file()
    call begin_operation(IO_OP_CLOSE)
    if (mpp_pe() == model_pes(1)) then
      call mpp_sync_self()
      call mpp_recv(acknowledgement, writer_pes(1), tag=protocol_tag)
      if (acknowledgement /= IO_OP_CLOSE) call mpp_error(FATAL, "offloading: invalid close acknowledgement")
    endif
    call end_operation()
    call mpp_sync()
    file_is_open = .false.
    session_finished = .true.
  end subroutine

  subroutine close_offloaded_files()
    call close_offloaded_file()
  end subroutine

  !> Enters main write loop for the offloaded writer PE
  !! to be called only by the writer after recv_domain_metadata.
  subroutine run_offload_writer()
    type(FmsNetcdfDomainFile_t) :: fileobj
    type(domain2d) :: domain
    integer :: operation, axis_size, naxes, direction, field_dims(2), axis_directions(2)
    character(len=string_length) :: name, cart, axis_names(2)
    real, allocatable :: data(:,:)

    call require_setup()
    if (mpp_pe() /= writer_pes(1)) call mpp_error(FATAL, "offloading: consumer must run on writer PE")
    if (.not. metadata_ready) call mpp_error(FATAL, "offloading: receive metadata before running writer")
    if (session_finished) call mpp_error(FATAL, "offloading: writer session already finished")
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    call mpp_define_domains((/1,global_size(1),1,global_size(2)/), (/1,1/), domain)
    call mpp_define_io_domain(domain, (/1,1/))
    allocate(data(global_size(1), global_size(2)))
    naxes = 0
    axis_directions = 0
    do
      call mpp_set_current_pelist(all_pes, no_sync=.true.)
      call mpp_recv(operation, model_pes(1), tag=protocol_tag)
      if (operation /= IO_OP_OPEN .and. .not. file_is_open) &
        call mpp_error(FATAL, "offloading: operation received before file open")
      select case (operation)
      case (IO_OP_OPEN)
        if (file_is_open) call mpp_error(FATAL, "offloading: file already open")
        call recv_string(name)
        call mpp_set_current_pelist(writer_pes, no_sync=.true.)
        if (.not. open_file(fileobj, trim(name), "overwrite", domain=domain)) &
          call mpp_error(FATAL, "offloading: unable to open "//trim(name))
        file_is_open = .true.
      case (IO_OP_REGISTER_AXIS)
        call recv_string(name)
        call recv_string(cart)
        call mpp_recv(axis_size, model_pes(1), tag=protocol_tag)
        if (naxes == 2) call mpp_error(FATAL, "offloading: only two axes are supported")
        select case (trim(cart))
        case ("x")
          direction = 1
        case ("y")
          direction = 2
        case default
          call mpp_error(FATAL, "offloading: invalid axis direction received")
        end select
        if (axis_size /= global_size(direction)) call mpp_error(FATAL, "offloading: invalid axis size received")
        if (any(axis_directions(:naxes) == direction)) call mpp_error(FATAL, "offloading: duplicate axis direction")
        naxes = naxes + 1
        axis_names(naxes) = name
        axis_directions(naxes) = direction
        call mpp_set_current_pelist(writer_pes, no_sync=.true.)
        call register_axis(fileobj, trim(name), trim(cart))
      case (IO_OP_REGISTER_FIELD)
        call recv_string(name)
        call mpp_recv(field_dims, 2, model_pes(1), tag=protocol_tag)
        if (any(field_dims < 1) .or. any(field_dims > naxes)) &
          call mpp_error(FATAL, "offloading: field refers to unregistered axes")
        if (any(axis_directions(field_dims) /= (/1,2/))) &
          call mpp_error(FATAL, "offloading: field dimensions must describe x then y")
        call mpp_set_current_pelist(writer_pes, no_sync=.true.)
        call register_field(fileobj, trim(name), "double", axis_names(field_dims))
      case (IO_OP_WRITE)
        call recv_string(name)
        call mpp_recv(data, size(data), model_pes(1), tag=protocol_tag)
        call mpp_set_current_pelist(writer_pes, no_sync=.true.)
        call write_data(fileobj, trim(name), data)
      case (IO_OP_GLOBAL_ATTRIBUTE)
        call recv_global_attribute(fileobj)
      case (IO_OP_CLOSE)
        call mpp_set_current_pelist(writer_pes, no_sync=.true.)
        call close_file(fileobj)
        file_is_open = .false.
        session_finished = .true.
        call mpp_set_current_pelist(all_pes, no_sync=.true.)
        call mpp_send(IO_OP_CLOSE, model_pes(1), tag=protocol_tag)
        call mpp_sync_self()
        exit
      case default
        call mpp_error(FATAL, "offloading: unknown operation received")
      end select
    enddo
    call mpp_set_current_pelist(writer_pes, no_sync=.true.)
    call mpp_deallocate_domain(domain)
  end subroutine

  !> Consumer half of global_metadata_offload. Called on the global pelist; returns on the writer pelist.
  subroutine recv_global_attribute(fileobj)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    character(len=string_length) :: name, str_value
    integer :: value_type, str_len
    real(kind=r4_kind) :: r4_value
    real(kind=r8_kind) :: r8_value
    integer(kind=i4_kind) :: i4_value
    integer(kind=i8_kind) :: i8_value

    call recv_string(name)
    call mpp_recv(value_type, model_pes(1), tag=protocol_tag)
    select case (value_type)
    case (ATTRIBUTE_R4)
      call mpp_recv(r4_value, model_pes(1), tag=protocol_tag)
      call mpp_set_current_pelist(writer_pes, no_sync=.true.)
      call register_global_attribute(fileobj, trim(name), r4_value)
    case (ATTRIBUTE_R8)
      call mpp_recv(r8_value, model_pes(1), tag=protocol_tag)
      call mpp_set_current_pelist(writer_pes, no_sync=.true.)
      call register_global_attribute(fileobj, trim(name), r8_value)
    case (ATTRIBUTE_I4)
      call mpp_recv(i4_value, model_pes(1), tag=protocol_tag)
      call mpp_set_current_pelist(writer_pes, no_sync=.true.)
      call register_global_attribute(fileobj, trim(name), i4_value)
    case (ATTRIBUTE_I8)
      call mpp_recv(i8_value, model_pes(1), tag=protocol_tag)
      call mpp_set_current_pelist(writer_pes, no_sync=.true.)
      call register_global_attribute(fileobj, trim(name), i8_value)
    case (ATTRIBUTE_STR)
      call mpp_recv(str_len, model_pes(1), tag=protocol_tag)
      call recv_string(str_value)
      call mpp_set_current_pelist(writer_pes, no_sync=.true.)
      call register_global_attribute(fileobj, trim(name), str_value, str_len=str_len)
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

  subroutine begin_operation(operation)
    integer, intent(in) :: operation

    call mpp_set_current_pelist(all_pes, no_sync=.true.)
    if (mpp_pe() == model_pes(1)) then
      call mpp_send(operation, writer_pes(1), tag=protocol_tag)
      call mpp_sync_self()
    endif
  end subroutine

  subroutine end_operation()
    ! mpp_send is nonblocking: finish sends before temporary buffers go away.
    if (mpp_pe() == model_pes(1)) call mpp_sync_self()
    call mpp_set_current_pelist(model_pes, no_sync=.true.)
  end subroutine

  subroutine check_string(value)
    character(len=*), intent(in) :: value
    if (len_trim(value) == 0 .or. len_trim(value) > string_length) &
      call mpp_error(FATAL, "offloading: names must contain 1 to 255 characters")
  end subroutine

  ! Encode characters using MPP integer transport, without assuming MPI_COMM_WORLD.
  subroutine send_string(value)
    character(len=*), intent(in) :: value
    integer :: buffer(string_length), i
    character(len=string_length) :: padded

    padded = value
    do i = 1, string_length
      buffer(i) = iachar(padded(i:i))
    enddo
    call mpp_send(buffer, string_length, writer_pes(1), tag=protocol_tag)
    call mpp_sync_self()
  end subroutine

  subroutine recv_string(value)
    character(len=string_length), intent(out) :: value
    integer :: buffer(string_length), i

    call mpp_recv(buffer, string_length, model_pes(1), tag=protocol_tag)
    do i = 1, string_length
      value(i:i) = achar(buffer(i))
    enddo
  end subroutine
#endif
end module offloading_io_mod
