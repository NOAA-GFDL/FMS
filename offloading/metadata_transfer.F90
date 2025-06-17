module metadata_transfer_mod
  use platform_mod
  use mpi
  use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_error, FATAL
  use fms_mod, only: string 

  implicit none
  public

  integer, parameter :: real8_type = 1
  integer, parameter :: real4_type = 2
  integer, parameter :: int8_type = 3
  integer, parameter :: int4_type = 4
  integer, parameter :: str_type = 5

  logical :: metadata_transfer_initialized = .false.

  integer :: metadata_transfer_type_mpi_id = -1 

  integer, parameter :: ATTR_NAME_MAX_LENGTH = 128 
  integer, parameter :: ATTR_VALUE_MAX_LENGTH = 128

  type :: metadata_type
    private
    integer                             :: attribute_type
    integer                             :: attribute_length
    character(len=ATTR_NAME_MAX_LENGTH) :: attribute_name
    real(kind=r8_kind)                  :: attribute_value(ATTR_VALUE_MAX_LENGTH)
    contains
      procedure :: fms_metadata_broadcast
      procedure :: get_attribute_type
      procedure :: set_attribute_type
      procedure :: get_attribute_value
      procedure :: set_attribute_value
      procedure :: get_attribute_name
      procedure :: set_attribute_name
  end type

contains

  !> Initialize the mpi datatype for future transfers
  subroutine fms_metadata_transfer_init()
    integer, dimension(0:3) :: lengths, types
    integer(KIND=MPI_ADDRESS_KIND), dimension(0:3) :: displacements
    integer :: ierror

    lengths = (/1, 1, ATTR_NAME_MAX_LENGTH, ATTR_VALUE_MAX_LENGTH/)
    displacements = (/0, sizeof(0), sizeof(0)*2, sizeof(0)*2 + sizeof(' ')*ATTR_NAME_MAX_LENGTH/) ! TODO, not sure if real size should be explicit
    types = (/MPI_INTEGER, MPI_INTEGER, MPI_CHARACTER, MPI_REAL/)
    call MPI_Type_create_struct(4, lengths, displacements, types, metadata_transfer_type_mpi_id, ierror)
    call MPI_Type_commit(metadata_transfer_type_mpi_id, ierror)
    if(ierror /= MPI_SUCCESS) then
      call mpp_error(FATAL, "fms_metadata_transfer_init: MPI_Type_create_struct failed")
    end if
    metadata_transfer_initialized = .true.
  end subroutine fms_metadata_transfer_init

  subroutine fms_metadata_broadcast(this)
    class(metadata_type), intent(inout) :: this
    integer :: ierror
    if (.not. metadata_transfer_initialized) then
      call mpp_error(FATAL, "fms_metadata_broadcast: metadata_transfer not initialized")
    end if

    ! Broadcast the metadata transfer type to all processes
    call MPI_Bcast(this, 1, metadata_transfer_type_mpi_id, mpp_root_pe(), MPI_COMM_WORLD, ierror)
    if (ierror /= MPI_SUCCESS) then
      call mpp_error(FATAL, "fms_metadata_broadcast: MPI_Bcast failed")
    end if

  end subroutine fms_metadata_broadcast

  subroutine fms_metadata_broadcast_all(metadata_objs)
    type(metadata_type), intent(inout) :: metadata_objs(:)
    integer :: ierror, i
    if (.not. metadata_transfer_initialized) then
      call mpp_error(FATAL, "fms_metadata_broadcast: metadata_transfer not initialized")
    end if

    do i=1, size(metadata_objs)
      call metadata_objs(i)%fms_metadata_broadcast()
    enddo

  end subroutine fms_metadata_broadcast_all

  ! Getter and Setter for attribute_type
  function get_attribute_type(this) result(val)
    class(metadata_type), intent(in) :: this
    integer :: val
    val = this%attribute_type
  end function

  subroutine set_attribute_type(this, val)
    class(metadata_type), intent(inout) :: this
    integer, intent(in) :: val
    this%attribute_type = val
  end subroutine

  ! Getter and Setter for attribute_value
  function get_attribute_value(this) result(val)
    class(metadata_type), intent(in) :: this
    real(kind=r8_kind) :: val(this%attribute_length)
    val = this%attribute_value(1:this%attribute_length)
  end function

  subroutine set_attribute_value(this, val)
    class(metadata_type), intent(inout) :: this
    real(kind=r8_kind), intent(in) :: val(:)
    if(size(val) .gt. ATTR_VALUE_MAX_LENGTH) then
      call mpp_error(FATAL, &
        "metadata_transfer_mod: attribute value array exceeds max length of "//string(ATTR_NAME_MAX_LENGTH))
    endif
    this%attribute_length = size(val)
    this%attribute_value(1:size(val)) = val
  end subroutine

  ! Getter and Setter for attribute_name
  function get_attribute_name(this) result(val)
    class(metadata_type), intent(in) :: this
    character(len=ATTR_NAME_MAX_LENGTH) :: val
    val = trim(this%attribute_name)
  end function

  subroutine set_attribute_name(this, val)
    class(metadata_type), intent(inout) :: this
    character(len=*), intent(in) :: val
    if(len(val) .gt. ATTR_NAME_MAX_LENGTH) then
      call mpp_error(FATAL, & 
        "metadata_transfer_mod: attribute name exceeds max length of "//string(ATTR_VALUE_MAX_LENGTH))
    endif
    this%attribute_name = val
  end subroutine

end module metadata_transfer_mod
