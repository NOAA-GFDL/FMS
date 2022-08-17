!***********************************************************************
!*           GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @author Ryan Mulhall
!> @email ryan.mulhall@noaa.gov
!! @brief Contains buffer types and routines for the diag manager
!!
!! @description Holds buffered data for fmsDiagVars_type objects
!! buffer0-5d types extend fmsDiagBuffer_class, and upon allocation
!! are added to the module's buffer_lists depending on it's dimension
!!
!! TODO: - the allocates assume 1 as the lower bound
!!       - add_to_buffer assumes copying to/from a 1:size(input array) 
module fms_diag_buffer_mod

use platform_mod
use iso_c_binding
use fms_diag_axis_object_mod, only: diagDomain_t
use time_manager_mod, only: time_type
use mpp_mod, only: mpp_error, FATAL
use diag_data_mod, only: DIAG_NULL, DIAG_NOT_REGISTERED

implicit none

private

!> @brief Object that holds buffered data and other diagnostics
!! Abstract to ensure use through its extensions(buffer0-5d types)
type, abstract :: fmsDiagBuffer_class
  class(*), dimension(:,:,:,:,:), allocatable :: remap_buffer !< remapped buffer data
  integer, allocatable, private               :: buffer_id !< index in buffer list
  !type(time_type), private                    :: init_time !< initialization time

  contains
  procedure :: get_remapped_buffer_pointer
  procedure :: flush_buffer
  procedure :: remap_buffer_data => remap_buffer
  procedure :: set_buffer_id
  ! TODO deferred routines, will require some interfaces 
  !procedure(allocate_buffer), deferred :: allocate_buffer
  !procedure, deferred :: get_buffer 
  !procedure, deferred :: initialize_buffer

end type fmsDiagBuffer_class

!> holds pointer to allocated buffer0-5d objects
type :: fmsDiagBufferContainer_type
  class(fmsDiagBuffer_class), allocatable :: diag_buffer_obj 
end type

!> Scalar buffer type to extend fmsDiagBufferContainer_type
type, extends(fmsDiagBuffer_class) :: buffer0d
  class(*), allocatable :: buffer !< scalar numberic buffer value
  contains
  procedure :: allocate_buffer => allocate_buffer_0d
  procedure :: get_buffer_data => get_buffer_0d
  procedure :: initialize_buffer => initialize_buffer_0d
  procedure :: add_to_buffer => add_to_buffer_0d

end type buffer0d

!> 1D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer1d
  class(*), allocatable :: buffer(:) !< 1D numeric data array
  contains
  procedure :: allocate_buffer => allocate_buffer_1d
  procedure :: get_buffer_data => get_buffer_1d
  procedure :: initialize_buffer => initialize_buffer_1d
  procedure :: add_to_buffer => add_to_buffer_1d
end type buffer1d

!> 2D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer2d
  class(*), allocatable :: buffer(:,:) !< 2D numeric data array
  contains
  procedure :: allocate_buffer => allocate_buffer_2d
  procedure :: get_buffer_data => get_buffer_2d
  procedure :: initialize_buffer => initialize_buffer_2d
  procedure :: add_to_buffer => add_to_buffer_2d
end type buffer2d

!> 3D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer3d
  class(*), allocatable :: buffer(:,:,:) !< 3D numeric data array
  contains
  procedure :: allocate_buffer => allocate_buffer_3d
  procedure :: get_buffer_data => get_buffer_3d
  procedure :: initialize_buffer => initialize_buffer_3d
  procedure :: add_to_buffer => add_to_buffer_3d
end type buffer3d

!> 4D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer4d
  class(*), allocatable :: buffer(:,:,:,:) !< 4D numeric data array
  contains
  procedure :: allocate_buffer => allocate_buffer_4d
  procedure :: get_buffer_data => get_buffer_4d
  procedure :: initialize_buffer => initialize_buffer_4d
  procedure :: add_to_buffer => add_to_buffer_4d
end type buffer4d

!> 5D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer5d
  class(*), allocatable :: buffer(:,:,:,:,:) !< 5D numeric data array
  contains
  procedure :: allocate_buffer => allocate_buffer_5d
  procedure :: get_buffer_data => get_buffer_5d
  procedure :: initialize_buffer => initialize_buffer_5d
  procedure :: add_to_buffer => add_to_buffer_5d
end type buffer5d

! public types
public :: buffer0d
public :: buffer1d
public :: buffer2d
public :: buffer3d
public :: buffer4d
public :: buffer5d
public :: fmsDiagBuffer_class
public :: fmsDiagBufferContainer_type
! public routines
public :: fms_diag_buffer_init 
! Module variables
logical,private :: module_is_initialized = .false. !< Flag indicating if the module is initialized
integer, parameter :: DEFAULT_BUFFER_LIST_SIZE = 64


logical, parameter, private :: DEBUG = .true. !< debugging output

! TODO interface with function declarations needed for deferred routines
!abstract interface
  !integer function allocate_buffer(buffobj, mold) &
  !result(rslt)
    !import buffer0d 
    !!import time_type 
    !class(fmsDiag), intent(inout) :: buffobj !< scalar buffer object
    !class(*),intent(in) :: mold !< allocates to the type of mold
    !!type(time_type), intent(in), optional :: init_time
  !end function
!end interface
!
contains


!!--------module routines

!> Initializes a list of diag buffers
!! and sets the buffer_id's to DIAG_NOT_REGISTERED
logical function fms_diag_buffer_init(buffobjs, buff_dims)
  type(fmsDiagBufferContainer_type), allocatable, intent(out) :: buffobjs(:) !< an array of buffer container types to allocate
  integer, intent(in)                                          :: buff_dims !< number of dimensions needed for the buffer data
  integer :: i

  allocate(buffobjs(DEFAULT_BUFFER_LIST_SIZE))
  do i=1, SIZE(buffobjs)
    buffobjs(i) = fms_diag_buffer_create_container(buff_dims)
    call buffobjs(i)%diag_buffer_obj%set_buffer_id(DIAG_NOT_REGISTERED)
  enddo
  fms_diag_buffer_init = allocated(buffobjs)
end function

!> creates a container type with a new (unallocated) buffer for the given dimensions
function fms_diag_buffer_create_container(buff_dims) &
result(rslt)
  integer, intent(in)                            :: buff_dims !< dimensions 
  type(fmsDiagBufferContainer_type), allocatable :: rslt

  allocate(rslt)
  select case (buff_dims)
    case (0)
      allocate(buffer0d :: rslt%diag_buffer_obj)
    case (1)
      allocate(buffer1d :: rslt%diag_buffer_obj)
    case (2)
      allocate(buffer2d :: rslt%diag_buffer_obj)
    case (3)
      allocate(buffer3d :: rslt%diag_buffer_obj)
    case (4)
      allocate(buffer4d :: rslt%diag_buffer_obj)
    case (5)
      allocate(buffer5d :: rslt%diag_buffer_obj)
    case default
      call mpp_error(FATAL, 'fms_diag_buffer_create_container: invalid number of dimensions given') 
  end select
end function 

!!--------generic routines for any fmsDiagBuffer_class objects 

!> setter for buffer_id for any buffer objects
subroutine set_buffer_id(this, id)
  class(fmsDiagBuffer_class), intent(inout) :: this !< buffer object to set id for
  integer, intent(in)                       :: id !< positive integer id to set
  if (.not.allocated(this%buffer_id) ) allocate(this%buffer_id)
  this%buffer_id = id
end subroutine

!> Remaps 0-5d data buffer from the given object onto a 5d array and sets it in the type
subroutine remap_buffer(this)
  class(fmsDiagBuffer_class), intent(inout) :: this !< any dimension buffer object
  integer                   :: buff_bounds(5)

  if( DEBUG) print *, 'remapping buffer'
  ! get num dimensions from type extension
  select type (this)
    type is (buffer0d)
      if (.not. allocated(this%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      if (DEBUG) print *, '0d buffer'
      ! get buffer data type to allocate and remap
      select type (buff => this%buffer)
        type is(integer(i4_kind))
          allocate (integer(kind=i4_kind) :: this%remap_buffer(1,1,1,1,1))
          this%remap_buffer = buff
        type is(integer(i8_kind))
          allocate (integer(kind=i8_kind) :: this%remap_buffer(1,1,1,1,1))
          this%remap_buffer = buff
        type is(real(r4_kind))
          allocate (real(kind=r4_kind) :: this%remap_buffer(1,1,1,1,1))
          this%remap_buffer = buff
        type is(real(r8_kind))
          allocate (real(kind=r8_kind) :: this%remap_buffer(1,1,1,1,1))
          this%remap_buffer = buff 
      end select
    type is (buffer1d)
      if (.not. allocated(this%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      if( DEBUG) print *, '1d buffer'
      buff_bounds(1) = SIZE(this%buffer, 1)
      select type (buff => this%buffer)
        type is(integer(i4_kind))
          allocate (integer(kind=i4_kind) :: this%remap_buffer(buff_bounds(1),1,1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
        type is(integer(i8_kind))
          allocate (integer(kind=i8_kind) :: this%remap_buffer(buff_bounds(1),1,1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
        type is(real(r4_kind))
          allocate (real(kind=r4_kind) :: this%remap_buffer(buff_bounds(1),1,1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
        type is(real(r8_kind))
          allocate (real(kind=r8_kind) :: this%remap_buffer(buff_bounds(1),1,1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
      end select
    type is (buffer2d)
      if (.not. allocated(this%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      buff_bounds(1) = SIZE(this%buffer, 1)
      buff_bounds(2) = SIZE(this%buffer, 2)
      select type (buff => this%buffer)
        type is(integer(i4_kind))
          allocate (integer(kind=i4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
        type is(integer(i8_kind))
          allocate (integer(kind=i8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
        type is(real(r4_kind))
          allocate (real(kind=r4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
        type is(real(r8_kind))
          allocate (real(kind=r8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
      end select
    type is (buffer3d)
      if (.not. allocated(this%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      buff_bounds(1) = SIZE(this%buffer, 1)
      buff_bounds(2) = SIZE(this%buffer, 2)
      buff_bounds(3) = SIZE(this%buffer, 3)
      select type (buff => this%buffer)
        type is(integer(i4_kind))
          allocate (integer(kind=i4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2), buff_bounds(3),1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
        type is(integer(i8_kind))
          allocate (integer(kind=i8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
        type is(real(r4_kind))
          allocate (real(kind=r4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
        type is(real(r8_kind))
          allocate (real(kind=r8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
      end select
    type is (buffer4d)
      buff_bounds(1) = SIZE(this%buffer, 1)
      buff_bounds(2) = SIZE(this%buffer, 2)
      buff_bounds(3) = SIZE(this%buffer, 3)
      buff_bounds(4) = SIZE(this%buffer, 4)
      select type (buff => this%buffer)
        type is(integer(i4_kind))
          allocate (integer(kind=i4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2), buff_bounds(3), buff_bounds(4),1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4),1 /))
        type is(integer(i8_kind))
          allocate (integer(kind=i8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2), buff_bounds(3), buff_bounds(4),1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4),1 /))
        type is(real(r4_kind))
          allocate (real(kind=r4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4),1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4),1 /))
        type is(real(r8_kind))
          allocate (real(kind=r8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4),1))
          this%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4),1 /))
      end select
    type is (buffer5d)
      buff_bounds(1) = SIZE(this%buffer, 1)
      buff_bounds(2) = SIZE(this%buffer, 2)
      buff_bounds(3) = SIZE(this%buffer, 3)
      buff_bounds(4) = SIZE(this%buffer, 4)
      buff_bounds(5) = SIZE(this%buffer, 5)
      select type (buff => this%buffer)
        type is(integer(i4_kind))
          allocate (integer(kind=i4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2), buff_bounds(3), buff_bounds(4), buff_bounds(5)))
          this%remap_buffer = buff 
        type is(integer(i8_kind))
          allocate (integer(kind=i8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2), buff_bounds(3), buff_bounds(4), buff_bounds(5)))
          this%remap_buffer = buff
        type is(real(r4_kind))
          allocate (real(kind=r4_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4), buff_bounds(5)))
          this%remap_buffer = buff
        type is(real(r8_kind))
          allocate (real(kind=r8_kind) :: this%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3), buff_bounds(4), buff_bounds(5)))
          this%remap_buffer = buff
      end select
    class default
      call mpp_error( FATAL, 'remap_buffer_pointer: invalid buffer type for remapping')
  end select
   
end subroutine 

!> @brief Gets the remapped buffer pointer
!! Will do the remapping if not set already
function get_remapped_buffer_pointer (this) &
result(rslt)
  class(fmsDiagBuffer_class), target, intent(inout) :: this !< any buffer object
  class(fmsDiagBuffer_class), pointer     :: objptr
  class(*), pointer, dimension(:,:,:,:,:)  :: rslt
  if(.not. allocated(this%remap_buffer)) then
    call remap_buffer(this) 
  endif
  objptr => this 
  rslt => objptr%remap_buffer
end function

!> Deallocates data fields from a buffer object
subroutine flush_buffer(this)
  class(fmsDiagBuffer_class), intent(inout) :: this !< any buffer object
  select type (this)
    type is (buffer0d)
      deallocate(this%buffer)
    type is (buffer1d)
      deallocate(this%buffer)
    type is (buffer2d)
      deallocate(this%buffer)
    type is (buffer3d)
      deallocate(this%buffer)
    type is (buffer4d)
      deallocate(this%buffer)
    type is (buffer5d)
      deallocate(this%buffer)
  end select
  deallocate(this%buffer_id)
end subroutine

!! -----------Type-specific routines for buffer0-5d 

!! allocations could be done in one routine for 0-5d if buffobj is changed to fmsDiagBuffer_class
!! not sure which approach would be better

!> allocates scalar buffer data to the given buff_type
subroutine allocate_buffer_0d(this, buff_type)
  class(buffer0d), intent(inout), target :: this !< scalar buffer object
  class(*),intent(in) :: buff_type !< allocates to the given type, value does not matter

  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer)
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer)
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer)
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer)
    class default
       call mpp_error("allocate_buffer_0d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

end subroutine

!> allocates 1D buffer data to given buff_type type
subroutine allocate_buffer_1d(this, buff_type, buff_size)
  class(buffer1d), intent(inout), target :: this !< scalar buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_size !< dimension bounds

  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_size))
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_size))
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_size))
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_size))
    class default
       call mpp_error("allocate_buffer_1d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

end subroutine
!> allocates a 2D buffer to given buff_type type
subroutine allocate_buffer_2d(this, buff_type, buff_sizes)
  class(buffer2d), intent(inout), target :: this !< 2D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(2) !< dimension sizes

  !if(.not. allocated(buffer_list_2d) .or. num_buffers(2)+1 .gt. BUFFER_LIST_SIZE) call allocate_buffer_list(2)
  !num_buffers(2) = num_buffers(2) + 1
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
    class default
       call mpp_error("allocate_buffer_1d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

end subroutine

!> allocates a 3D buffer to given buff_type type
subroutine allocate_buffer_3d(this, buff_type, buff_sizes)
  class(buffer3d), intent(inout), target :: this !< 3D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(3) !< dimension sizes 

  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
    class default
       call mpp_error("allocate_buffer_3d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select
end subroutine

!> allocates a 4D buffer to given buff_type type
subroutine allocate_buffer_4d(this, buff_type, buff_sizes)
  class(buffer4d), intent(inout), target :: this !< 4D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(4) !< dimension buff_sizes

  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
    class default
       call mpp_error("allocate_buffer_4d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select
end subroutine

!> allocates a 5D buffer to given buff_type type
subroutine allocate_buffer_5d(this, buff_type, buff_sizes)
  class(buffer5d), intent(inout), target :: this !< 5D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(5) !< dimension buff_sizes

  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),buff_sizes(5)))
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),buff_sizes(5)))
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),buff_sizes(5)))
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),buff_sizes(5)))
    class default
       call mpp_error("allocate_buffer_5d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select
end subroutine

!> @brief Gets buffer data from buffer0d type
!! @return copy of the buffer data
function get_buffer_0d (this) &
result(rslt)
  class (buffer0d), intent(in) :: this !< scalar buffer object
  class(*), allocatable :: rslt
  if (allocated(this%buffer)) then
    rslt = this%buffer
  else
    call mpp_error(FATAL, 'get_buffer_0d: buffer not allocated')
  endif
end function
!> @brief Gets buffer data from buffer1d type
!! @return copy of the buffer data
function get_buffer_1d (this) &
result(rslt)
  class (buffer1d), intent(in) :: this !< 1D buffer object
  class(*), allocatable :: rslt(:)
  if (allocated(this%buffer)) then
    rslt = this%buffer
  else
    call mpp_error(FATAL, 'get_buffer_1d: buffer not allocated')
  endif
end function
!> @brief Gets buffer data from buffer2d type
!! @return copy of the buffer data
function get_buffer_2d (this) &
result(rslt)
  class (buffer2d), intent(in) :: this !< 2D buffer object
  class(*), allocatable :: rslt(:,:)
  if (allocated(this%buffer)) then
    rslt = this%buffer
  else
    call mpp_error(FATAL, 'get_buffer_2d: buffer not allocated')
  endif
end function
!> @brief Gets buffer data from buffer3d type
!! @return copy of the buffer data
function get_buffer_3d (this) &
result(rslt)
  class (buffer3d), intent(in) :: this !< 3D buffer object
  class(*), allocatable :: rslt(:,:,:)
  if (allocated(this%buffer)) then
    rslt = this%buffer
  else
    call mpp_error(FATAL, 'get_buffer_3d: buffer not allocated')
  endif
end function
!> @brief Gets buffer data from buffer4d type
!! @return copy of the buffer data
function get_buffer_4d (this) &
result(rslt)
  class (buffer4d), intent(in) :: this !< 4D buffer object
  class(*), allocatable :: rslt(:,:,:,:)
  if (allocated(this%buffer)) then
    rslt = this%buffer
  else
    call mpp_error(FATAL, 'get_buffer_4d: buffer not allocated')
  endif
end function
!> @brief Gets buffer data from buffer5d type
!! @return copy of the buffer data
function get_buffer_5d (this) &
result(rslt)
  class (buffer5d), intent(in) :: this !< 5D buffer object
  class(*), allocatable :: rslt(:,:,:,:,:)
  if (allocated(this%buffer)) then
    rslt = this%buffer
  else
    call mpp_error(FATAL, 'get_buffer_5d: buffer not allocated')
  endif
end function

!> @brief Initializes a buffer to a given fill value
!! TODO default to 0?
subroutine initialize_buffer_0d (this, fillval)
  class(buffer0d), intent(inout) :: this !< scalar buffer object
  class(*) :: fillval

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_0d:' // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_0d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_1d (this, fillval)
  class(buffer1d), intent(inout) :: this !< 1D buffer object
  class(*) :: fillval !< fill value, must be same type as the allocated buffer in this

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_1d:' // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_1d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_2d (this, fillval)
  class(buffer2d), intent(inout) :: this !< 2D buffer object
  class(*) :: fillval

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_2d:' // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_2d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_3d (this, fillval)
  class(buffer3d), intent(inout) :: this !< 3D buffer object
  class(*) :: fillval

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_3d:' // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_3d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_4d (this, fillval)
  class(buffer4d), intent(inout) :: this !< allocated 4D buffer object
  class(*) :: fillval

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_4d:' // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_4d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_5d (this, fillval)
  class(buffer5d), intent(inout) :: this !< allocated 5D buffer object
  class(*) :: fillval

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_5d:' // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_5d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine

!> @brief Add values to 0d buffer
!! this will just call the init routine since there's only one value
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_0d(this, input_data)
  class(buffer0d), intent(inout) :: this !< allocated scalar buffer object
  class(*)      :: input_data
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_1d: buffer not yet allocated')
  call this%initialize_buffer(input_data)
end subroutine

!> @brief Add values to 1d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_1d(this, input_data)
  class(buffer1d), intent(inout) :: this !< 
  class(*)             :: input_data(:)
  integer            :: n !< number of elements in input data
  logical            :: type_error = .false.
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_1d: buffer not yet allocated')
  n = SIZE(input_data)
  if( n .gt. SIZE(this%buffer)) call mpp_error( FATAL,"add_to_buffer_1d: input data larger than allocated buffer")
  ! have to check both types for assignment
  select type( buffer => this%buffer )
  type is(integer(i4_kind))
    select type(input_data)
    type is(integer(i4_kind))
      buffer(1:n) = input_data(1:n)
    class default
      type_error = .true.
    end select
  type is(integer(i8_kind))
    select type(input_data)
    type is(integer(i8_kind))
      buffer(1:n) = input_data(1:n)
    class default
      type_error = .true.
    end select
  type is(real(r4_kind))
    select type(input_data)
    type is(real(r4_kind))
      buffer(1:n) = input_data(1:n)
    class default
      type_error = .true.
    end select
  type is(real(r8_kind))
    select type(input_data)
    type is(real(r8_kind))
      buffer(1:n) = input_data(1:n)
    class default
      type_error = .true.
    end select
  end select
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 2d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_2d(this, input_data)
  class(buffer2d), intent(inout) :: this
  class(*)             :: input_data(:,:)
  integer            :: n1, n2 !< number of elements per dimension
  logical            :: type_error
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_2d: buffer not yet allocated')
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2)) then
    call mpp_error( FATAL,"add_to_buffer_2d: input data larger than allocated buffer")
  endif
  !this%buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
  ! have to check both types for assignment
  select type( buffer => this%buffer )
  type is(integer(i4_kind))
    select type(input_data)
    type is(integer(i4_kind))
      buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
    class default
      type_error = .true.
    end select
  type is(integer(i8_kind))
    select type(input_data)
    type is(integer(i8_kind))
      buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
    class default
      type_error = .true.
    end select
  type is(real(r4_kind))
    select type(input_data)
    type is(real(r4_kind))
      buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
    class default
      type_error = .true.
    end select
  type is(real(r8_kind))
    select type(input_data)
    type is(real(r8_kind))
      buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
    class default
      type_error = .true.
    end select
  end select
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 3d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_3d(this, input_data)
  class(buffer3d), intent(inout) :: this
  class(*)             :: input_data(:,:,:)
  integer            :: n1, n2, n3 !< number of elements per dimension
  logical            :: type_error
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_3d: buffer not yet allocated')
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  n3 = SIZE(input_data, 3)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2) .or. &
    n3 .gt. SIZE(this%buffer, 3)) then
    call mpp_error( FATAL,"add_to_buffer_2d: input data larger than allocated buffer")
  endif
  !this%buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
  ! have to check both types for assignment
  select type( buffer => this%buffer )
  type is(integer(i4_kind))
    select type(input_data)
    type is(integer(i4_kind))
      buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
    class default
      type_error = .true.
    end select
  type is(integer(i8_kind))
    select type(input_data)
    type is(integer(i8_kind))
      buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
    class default
      type_error = .true.
    end select
  type is(real(r4_kind))
    select type(input_data)
    type is(real(r4_kind))
      buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
    class default
      type_error = .true.
    end select
  type is(real(r8_kind))
    select type(input_data)
    type is(real(r8_kind))
      buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
    class default
      type_error = .true.
    end select
  end select
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 4d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_4d(this, input_data)
  class(buffer4d), intent(inout) :: this
  class(*)             :: input_data(:,:,:,:)
  integer            :: n1, n2, n3, n4!< number of elements per dimension
  logical            :: type_error
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_4d: buffer not yet allocated')
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  n3 = SIZE(input_data, 3)
  n4 = SIZE(input_data, 4)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2) .or. &
    n3 .gt. SIZE(this%buffer, 3) .or. n4 .gt. SIZE(this%buffer, 4)) then
    call mpp_error( FATAL,"add_to_buffer_4d: input data larger than allocated buffer")
  endif
  ! have to check both types for assignment
  select type( buffer => this%buffer )
  type is(integer(i4_kind))
    select type(input_data)
    type is(integer(i4_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
    class default
      type_error = .true.
    end select
  type is(integer(i8_kind))
    select type(input_data)
    type is(integer(i8_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
    class default
      type_error = .true.
    end select
  type is(real(r4_kind))
    select type(input_data)
    type is(real(r4_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
    class default
      type_error = .true.
    end select
  type is(real(r8_kind))
    select type(input_data)
    type is(real(r8_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
    class default
      type_error = .true.
    end select
  end select
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_4d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 5d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_5d(this, input_data)
  class(buffer5d), intent(inout) :: this
  class(*)             :: input_data(:,:,:,:,:)
  integer            :: n1, n2, n3, n4, n5 !< number of elements per dimension
  logical            :: type_error
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_5d: buffer not yet allocated')
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  n3 = SIZE(input_data, 3)
  n4 = SIZE(input_data, 4)
  n5 = SIZE(input_data, 5)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2) .or. &
    n3 .gt. SIZE(this%buffer, 3) .or. n4 .gt. SIZE(this%buffer, 4) .or. &
    n5 .gt. SIZE(this%buffer, 5)) then
    call mpp_error( FATAL,"add_to_buffer_4d: input data larger than allocated buffer")
  endif
  !this%buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
  ! have to check both types for assignment
  select type( buffer => this%buffer )
  type is(integer(i4_kind))
    select type(input_data)
    type is(integer(i4_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
    class default
      type_error = .true.
    end select
  type is(integer(i8_kind))
    select type(input_data)
    type is(integer(i8_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
    class default
      type_error = .true.
    end select
  type is(real(r4_kind))
    select type(input_data)
    type is(real(r4_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
    class default
      type_error = .true.
    end select
  type is(real(r8_kind))
    select type(input_data)
    type is(real(r8_kind))
      buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
    class default
      type_error = .true.
    end select
  end select
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_5d: mismatch between allocated buffer and input data types')
end subroutine
end module fms_diag_buffer_mod
