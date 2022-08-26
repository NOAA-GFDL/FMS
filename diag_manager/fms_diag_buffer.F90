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
module fms_diag_buffer_mod

use platform_mod
use iso_c_binding
use fms_diag_axis_object_mod, only: diagDomain_t
use time_manager_mod, only: time_type
use mpp_mod, only: mpp_error, FATAL
use diag_data_mod, only: DIAG_NULL, DIAG_NOT_REGISTERED

implicit none

private

#ifdef use_yaml
!> @brief Object that holds buffered data and other diagnostics
!! Abstract to ensure use through its extensions(buffer0-5d types)
type, abstract :: fmsDiagBuffer_class
  integer, allocatable, private               :: buffer_id !< index in buffer list

  contains
  procedure :: flush_buffer
  procedure :: remap_buffer
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
  class(*), allocatable :: buffer(:) !< "scalar" numberic buffer value
                                     !! will only be allocated to hold 1 value
  character, private :: dims = '0'
  contains
  procedure :: allocate_buffer => allocate_buffer_0d
  procedure :: get_buffer_data => get_buffer_0d
  procedure :: initialize_buffer => initialize_buffer_0d
  procedure :: add_to_buffer => add_to_buffer_0d

end type buffer0d

!> 1D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer1d
  class(*), allocatable :: buffer(:) !< 1D numeric data array
  character, private :: dims = '1'
  contains
  procedure :: allocate_buffer => allocate_buffer_1d
  procedure :: get_buffer_data => get_buffer_1d
  procedure :: initialize_buffer => initialize_buffer_1d
  procedure :: add_to_buffer => add_to_buffer_1d
end type buffer1d

!> 2D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer2d
  class(*), allocatable :: buffer(:,:) !< 2D numeric data array
  character, private :: dims = '2'
  contains
  procedure :: allocate_buffer => allocate_buffer_2d
  procedure :: get_buffer_data => get_buffer_2d
  procedure :: initialize_buffer => initialize_buffer_2d
  procedure :: add_to_buffer => add_to_buffer_2d
end type buffer2d

!> 3D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer3d
  class(*), allocatable :: buffer(:,:,:) !< 3D numeric data array
  character, private :: dims = '3'
  contains
  procedure :: allocate_buffer => allocate_buffer_3d
  procedure :: get_buffer_data => get_buffer_3d
  procedure :: initialize_buffer => initialize_buffer_3d
  procedure :: add_to_buffer => add_to_buffer_3d
end type buffer3d

!> 4D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer4d
  class(*), allocatable :: buffer(:,:,:,:) !< 4D numeric data array
  character, private :: dims = '4'
  contains
  procedure :: allocate_buffer => allocate_buffer_4d
  procedure :: get_buffer_data => get_buffer_4d
  procedure :: initialize_buffer => initialize_buffer_4d
  procedure :: add_to_buffer => add_to_buffer_4d
end type buffer4d

!> 5D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer5d
  class(*), allocatable :: buffer(:,:,:,:,:) !< 5D numeric data array
  character, private :: dims = '5'
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
logical, parameter, private :: DEBUG = .true. !< debugging output

! TODO interface with function declarations needed for deferred routines
!abstract interface
  !integer function allocate_buffer(buffobj, buff_type) &
  !result(rslt)
    !import buffer0d
    !!import time_type
    !class(fmsDiag), intent(inout) :: buffobj !< scalar buffer object
    !class(*),intent(in) :: buff_type !< allocates to the type of buff_type
    !!type(time_type), intent(in), optional :: init_time
  !end function
!end interface
!
contains


!!--------module routines

!> Initializes a list of diag buffers
logical function fms_diag_buffer_init(buffobjs, buff_list_size)
  type(fmsDiagBufferContainer_type), allocatable, intent(out) :: buffobjs(:) !< an array of buffer container types to allocate
  integer, intent(in)                                          :: buff_list_size!< number of dimensions needed for the buffer data

  allocate(buffobjs(buff_list_size))
  fms_diag_buffer_init = allocated(buffobjs)
end function fms_diag_buffer_init

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
end function fms_diag_buffer_create_container

!!--------generic routines for any fmsDiagBuffer_class objects

!> setter for buffer_id for any buffer objects
subroutine set_buffer_id(this, id)
  class(fmsDiagBuffer_class), intent(inout) :: this !< buffer object to set id for
  integer, intent(in)                       :: id !< positive integer id to set
  if (.not.allocated(this%buffer_id) ) allocate(this%buffer_id)
  this%buffer_id = id
end subroutine set_buffer_id

!> Remaps 0-5d data buffer from the given object onto a 5d array pointer
function remap_buffer(buffobj)
  class(fmsDiagBuffer_class), target, intent(inout) :: buffobj !< any dimension buffer object
  class(*), pointer                                 :: remap_buffer(:,:,:,:,:)

  if( DEBUG) print *, 'remapping buffer'

  ! get num dimensions from type extension
  select type (buffobj)
    type is (buffer0d)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:1, 1:1, 1:1, 1:1) => buffobj%buffer
    type is (buffer1d)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:1, 1:1, 1:1, 1:1) => buffobj%buffer(1:size(buffobj%buffer,1))
    type is (buffer2d)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:1, 1:1, 1:1) => buffobj%buffer(:,:)
    type is (buffer3d)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), 1:1, 1:1) => buffobj%buffer(:,:,:)
    type is (buffer4d)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), &
                   1:size(buffobj%buffer,4), 1:1) => buffobj%buffer(:,:,:,:)
    type is (buffer5d)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), &
                   1:size(buffobj%buffer,4), 1:size(buffobj%buffer,5)) => buffobj%buffer(:,:,:,:,:)
    class default
      call mpp_error( FATAL, 'remap_buffer_pointer: invalid buffer type for remapping')
  end select

end function remap_buffer

!> Deallocates data fields from a buffer object
subroutine flush_buffer(this)
  class(fmsDiagBuffer_class), intent(inout) :: this !< any buffer object
  select type (this)
    type is (buffer0d)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer1d)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer2d)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer3d)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer4d)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer5d)
      if (allocated(this%buffer)) deallocate(this%buffer)
  end select
  if (allocated(this%buffer_id)) deallocate(this%buffer_id)
end subroutine flush_buffer

!! -----------Type-specific routines for buffer0-5d

!! allocations could be done in one routine for 0-5d if buffobj is changed to fmsDiagBuffer_class
!! not sure which approach would be better

!> allocates scalar buffer data to the given buff_type
subroutine allocate_buffer_0d(this, buff_type)
  class(buffer0d), intent(inout), target :: this !< scalar buffer object
  class(*),intent(in) :: buff_type !< allocates to the given type, value does not matter

  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(1))
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(1))
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(1))
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(1))
    class default
       call mpp_error("allocate_buffer_0d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

end subroutine allocate_buffer_0d

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

end subroutine allocate_buffer_1d
!> allocates a 2D buffer to given buff_type type
!! TODO fails with gnu
subroutine allocate_buffer_2d(this, buff_type, buff_sizes)
  class(buffer2d), intent(inout), target :: this !< 2D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(2) !< dimension sizes

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

end subroutine allocate_buffer_2d

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
end subroutine allocate_buffer_3d

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
end subroutine allocate_buffer_4d

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
end subroutine allocate_buffer_5d

!> @brief Gets buffer data from buffer0d type
!! @return copy of the buffer data
function get_buffer_0d (this) &
result(rslt)
  class (buffer0d), intent(in) :: this !< scalar buffer object
  class(*), allocatable :: rslt
  if (allocated(this%buffer)) then
    rslt = this%buffer(1)
  else
    call mpp_error(FATAL, 'get_buffer_0d: buffer not allocated')
  endif
end function get_buffer_0d
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
end function get_buffer_1d
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
end function get_buffer_2d
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
end function get_buffer_3d
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
end function get_buffer_4d
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
end function get_buffer_5d

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_0d (this, fillval)
  class(buffer0d), intent(inout) :: this !< scalar buffer object
  class(*), intent(in) :: fillval !< fill value, must be same type as the allocated buffer in this

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

end subroutine initialize_buffer_0d

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_1d (this, fillval)
  class(buffer1d), intent(inout) :: this !< 1D buffer object
  class(*), intent(in) :: fillval !< fill value, must be same type as the allocated buffer in this

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

end subroutine initialize_buffer_1d

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_2d (this, fillval)
  class(buffer2d), intent(inout) :: this !< 2D buffer object
  class(*), intent(in) :: fillval !< fill value, must be same type as the allocated buffer in this

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

end subroutine initialize_buffer_2d

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_3d (this, fillval)
  class(buffer3d), intent(inout) :: this !< 3D buffer object
  class(*), intent(in) :: fillval!< fill value, must be same type as the allocated buffer in this

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

end subroutine initialize_buffer_3d

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_4d (this, fillval)
  class(buffer4d), intent(inout) :: this !< allocated 4D buffer object
  class(*), intent(in) :: fillval!< fill value, must be same type as the allocated buffer in this

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

end subroutine initialize_buffer_4d

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_5d (this, fillval)
  class(buffer5d), intent(inout) :: this !< allocated 5D buffer object
  class(*), intent(in) :: fillval!< fill value, must be same type as the allocated buffer in this

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

end subroutine initialize_buffer_5d

!> @brief Add values to 0d buffer
!! this will just call the init routine since there's only one value
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_0d(this, input_data)
  class(buffer0d), intent(inout) :: this !< allocated scalar buffer object
  class(*), intent(in)      :: input_data !< data to copy into buffer
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_1d: buffer not yet allocated')
  call this%initialize_buffer(input_data)
end subroutine add_to_buffer_0d

!> @brief Copy values ( from 1 to size(input_data)) into a 1d buffer object
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_1d(this, input_data)
  class(buffer1d), intent(inout) :: this !< allocated 1d buffer object
  class(*), intent(in)            :: input_data(:) !< data to copy into the buffer
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
end subroutine add_to_buffer_1d

!> @brief Copy values ( from 1 to size(input_data)) into a 2d buffer object
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_2d(this, input_data)
  class(buffer2d), intent(inout) :: this !< allocated 2d buffer object
  class(*), intent(in)             :: input_data(:,:) !< 2d data array to copy into buffer
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
end subroutine add_to_buffer_2d

!> @brief Copy values ( from 1 to size(input_data)) into a 3d buffer object
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_3d(this, input_data)
  class(buffer3d), intent(inout) :: this !< allocated 3d buffer object
  class(*), intent(in)             :: input_data(:,:,:)!< 3d data array to copy into buffer
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
end subroutine add_to_buffer_3d

!> @brief Copy values ( from 1 to size(input_data)) into a 4d buffer object
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_4d(this, input_data)
  class(buffer4d), intent(inout) :: this !< allocated 4d buffer object
  class(*), intent(in)             :: input_data(:,:,:,:) !< 4d data to copy into buffer
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
end subroutine add_to_buffer_4d

!> @brief Copy values (from 1 to size(input_data)) into a 5d buffer object
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_5d(this, input_data)
  class(buffer5d), intent(inout) :: this !< allocated 5d buffer object
  class(*), intent(in)             :: input_data(:,:,:,:,:) !< 5d data to copy into buffer
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
end subroutine add_to_buffer_5d
#endif
end module fms_diag_buffer_mod
