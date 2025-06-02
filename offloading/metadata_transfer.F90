module metadata_transfer_mod
  use platform_mod
  use netcdf
  use mpi
  use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_error, FATAL

  implicit none
  public

  integer, parameter :: real8_type = 1
  integer, parameter :: real4_type = 2
  integer, parameter :: int8_type = 3
  integer, parameter :: int4_type = 4
  integer, parameter :: str_type = 5

  logical :: metadata_transfer_initialized = .false.

  integer :: metadata_transfer_type_mpi_id = -1 

  type :: metadata_type
    private
    integer            :: attribute_type
    real(kind=r8_kind) :: attribute_value(2)
    character(len=50)  :: attribute_name
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
    integer, dimension(0:2) :: lengths, types
    integer(KIND=MPI_ADDRESS_KIND), dimension(0:2) :: displacements
    integer :: ierror

    lengths = (/1, 2, 50/)
    displacements = (/0, sizeof(0), sizeof(0) + sizeof(0.0)/) ! TODO, not sure if real size should be explicit
    types = (/MPI_INTEGER, MPI_REAL, MPI_CHARACTER/)
    call MPI_Type_create_struct(3, lengths, displacements, types, metadata_transfer_type_mpi_id, ierror)
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
    real(kind=r8_kind) :: val(2)
    val = this%attribute_value
  end function

  subroutine set_attribute_value(this, val)
    class(metadata_type), intent(inout) :: this
    real(kind=r8_kind), intent(in) :: val(2)
    this%attribute_value = val
  end subroutine

  ! Getter and Setter for attribute_name
  function get_attribute_name(this) result(val)
    class(metadata_type), intent(in) :: this
    character(len=50) :: val
    val = this%attribute_name
  end function

  subroutine set_attribute_name(this, val)
    class(metadata_type), intent(inout) :: this
    character(len=*), intent(in) :: val
    this%attribute_name = val
  end subroutine

end module metadata_transfer_mod
