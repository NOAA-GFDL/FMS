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
module metadata_transfer_mod
  use platform_mod
#ifdef use_libMPI
  use mpi,     only: MPI_Type_create_struct, MPI_Type_commit, MPI_INTEGER, MPI_CHARACTER, &
                       MPI_DOUBLE, MPI_FLOAT, MPI_INT, MPI_LONG_INT, MPI_SUCCESS, MPI_ADDRESS_KIND
#endif
  use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_error, FATAL, mpp_get_current_pelist, mpp_npes
  use fms_mod, only: string

  implicit none

  public

#ifdef use_libMPI
  external MPI_Bcast
#endif

  integer, parameter :: real8_type = 1 !< enumeration for real(kind=8) data type
  integer, parameter :: real4_type = 2 !< enumeration for real(kind=4) data type
  integer, parameter :: int8_type = 3 !< enumeration for integer(kind=8) data type
  integer, parameter :: int4_type = 4 !< enumeration for integer(kind=4) data type
  integer, parameter :: str_type = 5 !< enumeration for string data type

  integer, parameter :: ATTR_NAME_MAX_LENGTH = 128
  integer, parameter :: ATTR_VALUE_MAX_LENGTH = 128

  !> Base class for broadcasting netcdf attribute data as a struct, holds the common fields
  !! and routines for initializing the mpi datatype so that children classes can
  !! be broadcasted.
  type, abstract :: metadata_class
    private
    integer                             :: mpi_type_id = -1 !< MPI datatype id corresponding to this data objects data
                                                            !! -1 if not set
    integer                             :: attribute_length = -1 !< length of the attribute value array, -1 if not set
    character(len=ATTR_NAME_MAX_LENGTH) :: attribute_name !< name of the attribute to write
    contains
      procedure :: fms_metadata_broadcast
      procedure :: fms_metadata_transfer_init
      procedure :: get_attribute_name
      procedure :: set_attribute_name
  end type


  !> Metadata class for real(kind=8) attribute values
  type, extends(metadata_class) :: metadata_r8_type
    real(r8_kind) :: attribute_value(ATTR_VALUE_MAX_LENGTH)
    contains
      procedure :: get_attribute_value => get_attribute_r8_value
      procedure :: set_attribute_value => set_attribute_r8_value
  end type metadata_r8_type

  !> Metadata class for real(kind=4) attribute values
  type, extends(metadata_class) :: metadata_r4_type
    real(r4_kind) :: attribute_value(ATTR_VALUE_MAX_LENGTH)
    contains
      procedure :: get_attribute_value => get_attribute_r4_value
      procedure :: set_attribute_value => set_attribute_r4_value
  end type metadata_r4_type

  !> Metadata class for integer(kind=8) attribute values
  type, extends(metadata_class) :: metadata_i8_type
    integer(i8_kind) :: attribute_value(ATTR_VALUE_MAX_LENGTH)
    contains
      procedure :: get_attribute_value => get_attribute_i8_value
      procedure :: set_attribute_value => set_attribute_i8_value
  end type metadata_i8_type

  !> Metadata class for integer(kind=4) attribute values
  type, extends(metadata_class) :: metadata_i4_type
    integer(i4_kind) :: attribute_value(ATTR_VALUE_MAX_LENGTH)
    contains
      procedure :: get_attribute_value => get_attribute_i4_value
      procedure :: set_attribute_value => set_attribute_i4_value
  end type metadata_i4_type

  !> Metadata class for string attribute values
  type, extends(metadata_class) :: metadata_str_type
    character(len=ATTR_VALUE_MAX_LENGTH) :: attribute_value
    contains
      procedure :: get_attribute_value => get_attribute_str_value
      procedure :: set_attribute_value => set_attribute_str_value
  end type metadata_str_type

  contains

  !> Initialize the mpi datatype for future broadcasts
  !! The metadata object's functions (not subroutines) are stored as fields in memory,
  !! so they need to be included in the MPI struct declaration.
  subroutine fms_metadata_transfer_init(this, dtype)
    class(metadata_class), intent(inout) :: this !<metadata object to initialize for mpi communication using the struct
    integer, intent(in) :: dtype !< data type and kind for the metadata's value
                                 !! must be real8_type, real4_type, int8_type, int4_type, or str_type
    integer, dimension(0:4) :: lengths, types
#ifdef use_libMPI
    integer(KIND=MPI_ADDRESS_KIND), dimension(0:4) :: displacements
    integer :: ierror, mpi_id

    !! since the actual data array is at the end of the struct, displacements are the same for all types
    !displacements = (/ 0, sizeof(0), sizeof(0)*2, sizeof(0)*3, &
    !                  sizeof(0)*3 + sizeof(' ')*ATTR_NAME_MAX_LENGTH, &
    !                  sizeof(0)*4 + sizeof(' ')*ATTR_NAME_MAX_LENGTH, &
    !                  sizeof(0)*4 :+ sizeof(' ')*ATTR_NAME_MAX_LENGTH + sizeof(' ') &
    !                /)
    displacements(0) = 0_MPI_ADDRESS_KIND ! id start address
    displacements(1) = displacements(0) + sizeof(0) ! attribute_length start address
    displacements(2) = displacements(1) + sizeof(0) ! attribute_name start adress
    displacements(3) = displacements(2) + sizeof(' ')*ATTR_NAME_MAX_LENGTH ! get_attribute_name() start address
    displacements(4) = displacements(3) + sizeof(' ')*ATTR_NAME_MAX_LENGTH ! attribute_value start address 

    select case(dtype)
    case(real8_type)
      types = (/MPI_INTEGER, MPI_INTEGER, MPI_CHARACTER, MPI_CHARACTER, MPI_DOUBLE/)
    case(real4_type)
      types = (/MPI_INTEGER, MPI_INTEGER, MPI_CHARACTER, MPI_CHARACTER, MPI_FLOAT/)
    case(int4_type)
      types = (/MPI_INTEGER, MPI_INTEGER, MPI_CHARACTER, MPI_CHARACTER, MPI_INT/)
    case(int8_type)
      types = (/MPI_INTEGER, MPI_INTEGER, MPI_CHARACTER, MPI_CHARACTER, MPI_LONG_INT/)
    case(str_type)
      types = (/MPI_INTEGER, MPI_INTEGER, MPI_CHARACTER, MPI_CHARACTER, MPI_CHARACTER/)
    case default
        call mpp_error(FATAL, "fms_metadata_transfer_init:: given dtype argument contains a unsupported type")
    end select

    !lengths = (/1, 1, 1, ATTR_NAME_MAX_LENGTH, ATTR_VALUE_MAX_LENGTH/)
    lengths = (/1, 1, ATTR_NAME_MAX_LENGTH, ATTR_NAME_MAX_LENGTH, ATTR_VALUE_MAX_LENGTH/)

    call MPI_Type_create_struct(4, lengths, displacements, types, mpi_id, ierror)
    if(ierror /= MPI_SUCCESS) then
      call mpp_error(FATAL, "fms_metadata_transfer_init: MPI_Type_create_struct failed")
    end if
    call MPI_Type_commit(mpi_id, ierror)
    if(ierror /= MPI_SUCCESS) then
      call mpp_error(FATAL, "fms_metadata_transfer_init: MPI_Type_commit failed")
    end if
    this%mpi_type_id = mpi_id
#else
  call mpp_error(FATAL, "fms_metadata_transfer_init: MPI library not enabled, cannot initialize metadata transfer")
#endif
  end subroutine fms_metadata_transfer_init

  !> Broadcast the entire metadata object to all PEs in the current pelist
  subroutine fms_metadata_broadcast(this)
    class(metadata_class), intent(inout) :: this !< object that inherits metadata_class
    integer :: ierror, curr_comm_id
    integer, allocatable :: broadcasting_pes(:)
    if (this%mpi_type_id .eq. -1) then
      call mpp_error(FATAL, "fms_metadata_broadcast: metadata_transfer not initialized")
    end if

    allocate(broadcasting_pes(mpp_npes()))
    call mpp_get_current_pelist(broadcasting_pes, commID=curr_comm_id)

#ifdef use_libMPI
    ! Broadcast the metadata transfer type to all processes
    select type(this)
    type is (metadata_r8_type)
      call MPI_Bcast(this, 1, this%mpi_type_id, mpp_root_pe(), curr_comm_id, ierror)
    type is (metadata_r4_type)
      call MPI_Bcast(this, 1, this%mpi_type_id, mpp_root_pe(), curr_comm_id, ierror)
    type is (metadata_i4_type)
      call MPI_Bcast(this, 1, this%mpi_type_id, mpp_root_pe(), curr_comm_id, ierror)
    type is (metadata_i8_type)
      call MPI_Bcast(this, 1, this%mpi_type_id, mpp_root_pe(), curr_comm_id, ierror)
    type is (metadata_str_type)
      call MPI_Bcast(this, 1, this%mpi_type_id, mpp_root_pe(), curr_comm_id, ierror)
    end select
    if (ierror /= MPI_SUCCESS) then
      call mpp_error(FATAL, "fms_metadata_broadcast: MPI_Bcast failed")
    end if
#else
  call mpp_error(FATAL, "fms_metadata_broadcast: MPI library not enabled")
#endif


  end subroutine fms_metadata_broadcast

  !> Broadcast an array of metadata objects to all PEs in the current pelist
  subroutine fms_metadata_broadcast_all(metadata_objs)
    class(metadata_class), intent(inout) :: metadata_objs(:) !< list of metadata objects
    integer :: i

    do i=1, size(metadata_objs)
      if (metadata_objs(i)%mpi_type_id .eq. -1) then
        call mpp_error(FATAL, "fms_metadata_broadcast_all: metadata_transfer not initialized")
      end if
      call metadata_objs(i)%fms_metadata_broadcast()
    enddo

  end subroutine fms_metadata_broadcast_all

  !> Getter for real 8 attribute_value
  function get_attribute_r8_value(this) result(val)
    class(metadata_r8_type), intent(inout) :: this
    real(r8_kind), allocatable :: val(:)
    val = this%attribute_value(1:this%attribute_length)
  end function

  !> Setter for real 8 attribute_value
  subroutine set_attribute_r8_value(this, val)
    class(metadata_r8_type), intent(inout) :: this
    real(r8_kind), intent(in) :: val(:) !< 8 byte real value to set attribute value to
    if(size(val) .gt. ATTR_VALUE_MAX_LENGTH) then
      call mpp_error(FATAL, &
        "metadata_transfer_mod: attribute value array exceeds max length of "//string(ATTR_NAME_MAX_LENGTH))
    endif
    this%attribute_length = size(val)
    this%attribute_value(1:size(val)) = val
  end subroutine

  !> Getter for real 4 attribute_value
  function get_attribute_r4_value(this) result(val)
    class(metadata_r4_type), intent(inout) :: this
    real(r4_kind), allocatable :: val(:)
    val = this%attribute_value(1:this%attribute_length)
  end function

  !> Setter for real 4 attribute_value
  subroutine set_attribute_r4_value(this, val)
    class(metadata_r4_type), intent(inout) :: this
    real(r4_kind), intent(in) :: val(:) !< 4 byte real attribute to set
    if(size(val) .gt. ATTR_VALUE_MAX_LENGTH) then
      call mpp_error(FATAL, &
        "metadata_transfer_mod: attribute value array exceeds max length of "//string(ATTR_NAME_MAX_LENGTH))
    endif
    this%attribute_length = size(val)
    this%attribute_value(1:size(val)) = val
  end subroutine

  !> Getter for integer(kind=8) attribute_value
  function get_attribute_i8_value(this) result(val)
    class(metadata_i8_type), intent(inout) :: this
    integer(i8_kind), allocatable :: val(:)
    val = this%attribute_value(1:this%attribute_length)
  end function

  !> Setter for integer(kind=8) attribute_value
  subroutine set_attribute_i8_value(this, val)
    class(metadata_i8_type), intent(inout) :: this
    integer(i8_kind), intent(in) :: val(:) !< 8 byte int attribute to set
    if(size(val) .gt. ATTR_VALUE_MAX_LENGTH) then
      call mpp_error(FATAL, &
        "metadata_transfer_mod: attribute value array exceeds max length of "//string(ATTR_NAME_MAX_LENGTH))
    endif
    this%attribute_length = size(val)
    this%attribute_value(1:size(val)) = val
  end subroutine

  !> Getter for integer(kind=4) attribute_value
  function get_attribute_i4_value(this) result(val)
    class(metadata_i4_type), intent(inout) :: this
    integer(i4_kind), allocatable :: val(:)
    val = this%attribute_value(1:this%attribute_length)
  end function

  !> Setter for integer(kind=4) attribute_value
  subroutine set_attribute_i4_value(this, val)
    class(metadata_i4_type), intent(inout) :: this
    integer(i4_kind), intent(in) :: val(:) !< 4 byte integer to set attribute value to
    if(size(val) .gt. ATTR_VALUE_MAX_LENGTH) then
      call mpp_error(FATAL, &
        "metadata_transfer_mod: attribute value array exceeds max length of "//string(ATTR_NAME_MAX_LENGTH))
    endif
    this%attribute_length = size(val)
    this%attribute_value(1:size(val)) = val
  end subroutine

  !> Getter for string attribute_value
  function get_attribute_str_value(this) result(val)
    class(metadata_str_type), intent(inout) :: this
    character(len=:), allocatable :: val
    val = this%attribute_value(1:this%attribute_length)
  end function

  !> Setter for string attribute_value
  subroutine set_attribute_str_value(this, val)
    class(metadata_str_type), intent(inout) :: this
    character(len=*), intent(in) :: val !< character string to set attribute value to
    if(len(val) .gt. ATTR_VALUE_MAX_LENGTH) then
      call mpp_error(FATAL, &
        "metadata_transfer_mod: attribute value array exceeds max length of "//string(ATTR_NAME_MAX_LENGTH))
    endif
    this%attribute_length = len(val)
    this%attribute_value(1:len(val)) = val
  end subroutine

  !> Getter for attribute_name (for all metadata types)
  function get_attribute_name(this) result(val)
    class(metadata_class), intent(inout) :: this
    character(len=ATTR_NAME_MAX_LENGTH) :: val     val = trim(this%attribute_name)
  end function

  !> Setter for attribute_name (for all metadata types)
  subroutine set_attribute_name(this, val)
    class(metadata_class), intent(inout) :: this
    character(len=*), intent(in) :: val !< character string to set attribute name to
    if(len(val) .gt. ATTR_NAME_MAX_LENGTH) then
      call mpp_error(FATAL, &
        "metadata_transfer_mod: attribute name exceeds max length of "//string(ATTR_VALUE_MAX_LENGTH))
    endif
    this%attribute_name = val
  end subroutine

end module metadata_transfer_mod
