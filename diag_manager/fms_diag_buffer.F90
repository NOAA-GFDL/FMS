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
  integer, allocatable, public :: num_elements(:) !< used in time-averaging
  class(*), allocatable, public :: count_0d(:) !< used in time-averaging along with
                                       !! counter which is stored in the child types (bufferNd)
  character(len=2), public :: typestr !<set to allocated data type & kind value, one of i4,i8,r4,r8
  integer, allocatable, public :: buffer_dims(:) !< holds the size of each dimension in the buffer
  contains

  procedure :: flush_buffer
  procedure :: remap_buffer
  procedure :: set_buffer_id
  ! TODO deferred routines, will require some interfaces
  !procedure(allocate_buffer), deferred :: allocate_buffer
  !procedure, deferred :: get_buffer
  !procedure, deferred :: initialize_buffer

end type fmsDiagBuffer_class

!> holds an allocated buffer0-5d object
type :: fmsDiagBufferContainer_type
  class(fmsDiagBuffer_class), allocatable :: diag_buffer_obj !< any 0-5d buffer object
end type

!> Scalar buffer type to extend fmsDiagBufferContainer_type
type, extends(fmsDiagBuffer_class) :: buffer0d_type
  class(*), allocatable :: buffer(:) !< "scalar" numberic buffer value
                                     !! will only be allocated to hold 1 value
  class(*), allocatable :: counter(:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_0d
  procedure :: initialize_buffer => initialize_buffer_0d
  procedure :: add_to_buffer => add_to_buffer_0d
  generic :: get_buffer => get_0d_int4, get_0d_int8, get_0d_real4, get_0d_real8
  procedure, private :: get_0d_int4, get_0d_int8, get_0d_real4, get_0d_real8

end type buffer0d_type

!> 1D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer1d_type
  class(*), allocatable :: buffer(:) !< 1D numeric data array
  class(*), allocatable :: counter(:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_1d
  procedure :: initialize_buffer => initialize_buffer_1d
  procedure :: add_to_buffer => add_to_buffer_1d
  generic :: get_buffer => get_1d_int4, get_1d_int8, get_1d_real4, get_1d_real8
  procedure, private :: get_1d_int4, get_1d_int8, get_1d_real4, get_1d_real8
end type buffer1d_type

!> 2D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer2d_type
  class(*), allocatable :: buffer(:,:) !< 2D numeric data array
  class(*), allocatable :: counter(:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_2d
  procedure :: initialize_buffer => initialize_buffer_2d
  procedure :: add_to_buffer => add_to_buffer_2d
  generic :: get_buffer => get_2d_int4, get_2d_int8, get_2d_real4, get_2d_real8
  procedure, private :: get_2d_int4, get_2d_int8, get_2d_real4, get_2d_real8
end type buffer2d_type

!> 3D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer3d_type
  class(*), allocatable :: buffer(:,:,:) !< 3D numeric data array
  class(*), allocatable :: counter(:,:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_3d
  procedure :: initialize_buffer => initialize_buffer_3d
  procedure :: add_to_buffer => add_to_buffer_3d
  generic :: get_buffer => get_3d_int4, get_3d_int8, get_3d_real4, get_3d_real8
  procedure, private :: get_3d_int4, get_3d_int8, get_3d_real4, get_3d_real8
end type buffer3d_type

!> 4D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer4d_type
  class(*), allocatable :: buffer(:,:,:,:) !< 4D numeric data array
  class(*), allocatable :: counter(:,:,:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_4d
  procedure :: initialize_buffer => initialize_buffer_4d
  procedure :: add_to_buffer => add_to_buffer_4d
  generic :: get_buffer => get_4d_int4, get_4d_int8, get_4d_real4, get_4d_real8
  procedure, private :: get_4d_int4, get_4d_int8, get_4d_real4, get_4d_real8
end type buffer4d_type

!> 5D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer5d_type
  class(*), allocatable :: buffer(:,:,:,:,:) !< 5D numeric data array
  class(*), allocatable :: counter(:,:,:,:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_5d
  procedure :: initialize_buffer => initialize_buffer_5d
  procedure :: add_to_buffer => add_to_buffer_5d
  generic :: get_buffer => get_5d_int4, get_5d_int8, get_5d_real4, get_5d_real8
  procedure, private :: get_5d_int4, get_5d_int8, get_5d_real4, get_5d_real8
end type buffer5d_type

! public types
public :: buffer0d_type
public :: buffer1d_type
public :: buffer2d_type
public :: buffer3d_type
public :: buffer4d_type
public :: buffer5d_type
public :: fmsDiagBuffer_class
public :: fmsDiagBufferContainer_type

! public routines
public :: fms_diag_buffer_init

contains

!!--------module routines

!> Initializes a list of diag buffers
logical function fms_diag_buffer_init(buffobjs, buff_list_size)
  type(fmsDiagBufferContainer_type), allocatable, intent(out) :: buffobjs(:) !< an array of buffer container types
                                                                             !! to allocate
  integer, intent(in)                                          :: buff_list_size !< number of dimensions needed for
                                                                             !! the buffer data

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
      allocate(buffer0d_type :: rslt%diag_buffer_obj)
    case (1)
      allocate(buffer1d_type :: rslt%diag_buffer_obj)
    case (2)
      allocate(buffer2d_type :: rslt%diag_buffer_obj)
    case (3)
      allocate(buffer3d_type :: rslt%diag_buffer_obj)
    case (4)
      allocate(buffer4d_type :: rslt%diag_buffer_obj)
    case (5)
      allocate(buffer5d_type :: rslt%diag_buffer_obj)
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

  ! get num dimensions from type extension
  select type (buffobj)
    type is (buffer0d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:1, 1:1, 1:1, 1:1) => buffobj%buffer
    type is (buffer1d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:1, 1:1, 1:1, 1:1) => buffobj%buffer(1:size(buffobj%buffer,1))
    type is (buffer2d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:1, 1:1, 1:1) => buffobj%buffer(:,:)
    type is (buffer3d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), 1:1, 1:1) => &
                                                                                          & buffobj%buffer(:,:,:)
    type is (buffer4d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), &
                   1:size(buffobj%buffer,4), 1:1) => buffobj%buffer(:,:,:,:)
    type is (buffer5d_type)
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
    type is (buffer0d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer1d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer2d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer3d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer4d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
    type is (buffer5d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
  end select
  if (allocated(this%buffer_id)) deallocate(this%buffer_id)
end subroutine flush_buffer

!! -----------Type-specific routines for buffer0-5d

!! allocations could be done in one routine for 0-5d if buffobj is changed to fmsDiagBuffer_class
!! not sure which approach would be better

!> allocates scalar buffer data to the given buff_type
subroutine allocate_buffer_0d(this, buff_type, diurnal_samples)
  class(buffer0d_type), intent(inout), target :: this !< scalar buffer object
  class(*),intent(in) :: buff_type !< allocates to the given type, value does not matter
  integer, optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_0d: buffer already allocated")
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(1))
      allocate(integer(kind=i4_kind) :: this%counter(1))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i4'
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(1))
      allocate(integer(kind=i8_kind) :: this%counter(1))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i8'
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(1))
      allocate(real(kind=r4_kind) :: this%counter(1))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r4'
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(1))
      allocate(real(kind=r8_kind) :: this%counter(1))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r8'
    class default
       call mpp_error("allocate_buffer_0d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

  allocate(this%num_elements(n_samples))
  allocate(this%buffer_dims(1))
  this%num_elements = 0
  this%count_0d   = 0
  this%buffer_dims(1) = 1

end subroutine allocate_buffer_0d

!> allocates 1D buffer data to given buff_type type
subroutine allocate_buffer_1d(this, buff_type, buff_size, diurnal_samples)
  class(buffer1d_type), intent(inout), target :: this !< scalar buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_size !< dimension bounds
  integer, optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_1d: buffer already allocated")
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_size))
      allocate(integer(kind=i4_kind) :: this%counter(buff_size))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i4'
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_size))
      allocate(integer(kind=i8_kind) :: this%counter(buff_size))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i8'
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_size))
      allocate(real(kind=r4_kind) :: this%count_0d(buff_size))
      allocate(real(kind=r4_kind) :: this%counter(n_samples))
      this%counter = 0
      this%typestr = 'r4'
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_size))
      allocate(real(kind=r8_kind) :: this%count_0d(buff_size))
      allocate(real(kind=r8_kind) :: this%counter(n_samples))
      this%counter = 0
      this%typestr = 'r8'
    class default
       call mpp_error("allocate_buffer_1d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

  allocate(this%num_elements(n_samples))
  allocate(this%buffer_dims(1))
  this%num_elements = 0
  this%count_0d   = 0
  this%buffer_dims(1) = buff_size

end subroutine allocate_buffer_1d
!> allocates a 2D buffer to given buff_type type
!! TODO fails with gnu
subroutine allocate_buffer_2d(this, buff_type, buff_sizes, diurnal_samples)
  class(buffer2d_type), intent(inout), target :: this !< 2D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(2) !< dimension sizes
  integer, optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_2d: buffer already allocated")
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i4'
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i8'
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r4'
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r8_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%typestr = 'r4'
    class default
       call mpp_error("allocate_buffer_1d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select
  allocate(this%num_elements(n_samples))
  allocate(this%buffer_dims(2))
  this%num_elements = 0
  this%count_0d   = 0
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)

end subroutine allocate_buffer_2d

!> allocates a 3D buffer to given buff_type type
subroutine allocate_buffer_3d(this, buff_type, buff_sizes, diurnal_samples)
  class(buffer3d_type), intent(inout), target :: this !< 3D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(3) !< dimension sizes
  integer, optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_3d: buffer already allocated")
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i4'
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i8'
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r4'
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r8_kind) :: this%counter( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%typestr = 'r4'
      this%counter = 0
    class default
       call mpp_error("allocate_buffer_3d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

  allocate(this%num_elements(n_samples))
  this%num_elements = 0
  this%count_0d   = 0
  allocate(this%buffer_dims(3))
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)
  this%buffer_dims(3) = buff_sizes(3)

end subroutine allocate_buffer_3d

!> allocates a 4D buffer to given buff_type type
subroutine allocate_buffer_4d(this, buff_type, buff_sizes, diurnal_samples)
  class(buffer4d_type), intent(inout), target :: this !< 4D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(4) !< dimension buff_sizes
  integer, optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_4d: buffer already allocated")
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i4'
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i8'
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r4'
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r8'
    class default
       call mpp_error("allocate_buffer_4d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select

  allocate(this%num_elements(n_samples))
  this%num_elements = 0
  this%count_0d   = 0
  allocate(this%buffer_dims(4))
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)
  this%buffer_dims(3) = buff_sizes(3)
  this%buffer_dims(4) = buff_sizes(4)

end subroutine allocate_buffer_4d

!> allocates a 5D buffer to given buff_type type
subroutine allocate_buffer_5d(this, buff_type, buff_sizes, diurnal_samples)
  class(buffer5d_type), intent(inout), target :: this !< 5D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(5) !< dimension buff_sizes
  integer, optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_5d: buffer already allocated")
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                                  & buff_sizes(5)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                   & buff_sizes(5)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i4'
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                  & buff_sizes(5)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                  & buff_sizes(5)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'i8'
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                               & buff_sizes(5)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                & buff_sizes(5)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r4'
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                               & buff_sizes(5)))
      allocate(real(kind=r8_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                & buff_sizes(5)))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%counter = 0
      this%typestr = 'r8'
    class default
       call mpp_error("allocate_buffer_5d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4",&
           FATAL)
  end select
  allocate(this%num_elements(n_samples))
  this%num_elements = 0
  this%count_0d   = 0
  allocate(this%buffer_dims(5))
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)
  this%buffer_dims(3) = buff_sizes(3)
  this%buffer_dims(4) = buff_sizes(4)
  this%buffer_dims(5) = buff_sizes(5)
end subroutine allocate_buffer_5d

!! gonna leave these for when we stop caring about gnu <11
!> @brief Gets buffer data from buffer0d_type type
!! @return copy of the buffer data
!function get_buffer_0d (this) &
!result(rslt)
  !class (buffer0d_type), intent(in) :: this !< scalar buffer object
  !class(*), allocatable :: rslt
  !if (allocated(this%buffer)) then
    !rslt = this%buffer(1)
  !else
    !call mpp_error(FATAL, 'get_buffer_0d: buffer not allocated')
  !endif
!end function get_buffer_0d

!> Gets real(r4_kind) buffer data from a scalar buffer object
!! called through buffobj%get_buffer(r4_outputdata)
subroutine get_0d_real4 (this, buff_out)
  class(buffer0d_type), intent(in) :: this !< 0d allocated buffer object
  real(r4_kind), allocatable, intent(out)  :: buff_out !< output of copied buffer data
  allocate(buff_out)
  select type (buff=>this%buffer)
    type is (real(r4_kind))
      buff_out = buff(1)
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!> Gets real(r8_kind) buffer data from a 0d buffer object
!! called through this%get_buffer(r8_outputdata)
subroutine get_0d_real8 (this, buff_out)
  class(buffer0d_type), intent(in) :: this !< 0d allocated buffer objects
  real(r8_kind), allocatable, intent(out)  :: buff_out !< output of copied buffer data
  allocate(buff_out)
  select type (buff=>this%buffer)
    type is (real(r8_kind))
      buff_out = buff(1)
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!> Gets integer(i4_kind) buffer data from 0d buffer object
!> called through this%get_buffer(i4_outputdata)
subroutine get_0d_int4 (this, buff_out)
  class(buffer0d_type), intent(in) :: this !< 0d allocated buffer objects
  integer(i4_kind), allocatable, intent(out)  :: buff_out !< output of copied buffer data
  allocate(buff_out)
  select type (buff=>this%buffer)
    type is (integer(r4_kind))
      buff_out = buff(1)
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!> Gets integer(i8_kind) buffer data from 0d buffer object
!> called through this%get_buffer(i8_outputdata)
subroutine get_0d_int8 (this, buff_out)
  class(buffer0d_type), intent(in) :: this !< 0d allocated buffer objects
  integer(i8_kind), allocatable, intent(out)  :: buff_out !< output of copied buffer data
  allocate(buff_out)
  select type (buff=>this%buffer)
    type is (integer(r8_kind))
      buff_out = buff(1)
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!> @brief Gets buffer data from buffer1d_type type
!! @return copy of the buffer data
!function get_buffer_1d (this) &
!result(rslt)
  !class (buffer1d_type), target, intent(in) :: this !< 1D buffer object
  !class(*), allocatable :: rslt(:)
  !select type(buff=>this%buffer)
    !type is(integer(8))
      !allocate(integer(8) :: rslt(size(buff)))
    !type is(integer(4))
      !allocate(integer(4) :: rslt(size(buff)))
    !type is(real(4))
      !allocate(real(4) :: rslt(size(buff)))
    !type is(real(8))
      !allocate(real(8) :: rslt(size(buff)))
  !end select
  !if (allocated(this%buffer)) then
    !rslt = this%buffer
  !else
    !call mpp_error(FATAL, 'get_buffer_1d: buffer not allocated')
  !endif
!end function get_buffer_1d

!> Gets real(r4_kind) buffer data from a 1d buffer object
!! called through this%get_buffer(r4_outputdata)
subroutine get_1d_real4 (this, buff_out)
  class(buffer1d_type), intent(in) :: this !< 1d allocated buffer object
  real(r4_kind), allocatable, intent(out)  :: buff_out(:) !< output of copied buffer data
                                            !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer)))
  select type (buff=>this%buffer)
    type is (real(r4_kind))
      buff_out = buff
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!> Gets real(r8_kind) buffer data from a 1d buffer object
!! called through this%get_buffer(r8_outputdata)
subroutine get_1d_real8 (this, buff_out)
  class(buffer1d_type), intent(in) :: this !< 1d allocated buffer object
  real(r8_kind), allocatable, intent(out)  :: buff_out(:) !< output of copied buffer data
                                         !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer)))
  select type (buff=>this%buffer)
    type is (real(r8_kind))
      buff_out = buff
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!> Gets integer(i4_kind) buffer data from 1d buffer object
!> called through this%get_buffer(i4_outputdata)
subroutine get_1d_int4 (this, buff_out)
  class(buffer1d_type), intent(in) :: this !< 1d allocated buffer object
  integer(i4_kind), allocatable, intent(out)  :: buff_out(:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer)))
  select type (buff=>this%buffer)
    type is (integer(r4_kind))
      buff_out = buff
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!> Gets integer(i8_kind) buffer data from 1d buffer object
!> called through this%get_buffer(i8_outputdata)
subroutine get_1d_int8 (this, buff_out)
  class(buffer1d_type), intent(in) :: this !< 1d allocated buffer object
  integer(i8_kind), allocatable, intent(out)  :: buff_out(:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer)))
  select type (buff=>this%buffer)
    type is (integer(r8_kind))
      buff_out = buff
    class default
      call mpp_error(FATAL, "incorrect output type for allocated buffer")
  end select
end subroutine

!! breaks gcc 10 and prior
!function get_buffer_2d (this) &
!result(rslt)
  !class (buffer2d_type), intent(in) :: this !< 2D buffer object
  !class(*), allocatable :: rslt(:,:)
  !select type(buff=>this%buffer)
    !type is(integer(8))
      !allocate(integer(8) :: rslt(size(buff,1), size(buff,2)))
    !type is(integer(4))
      !allocate(integer(4) :: rslt(size(buff,1), size(buff,2)))
    !type is(real(4))
      !allocate(real(4) :: rslt(size(buff,1), size(buff,2)))
    !type is(real(8))
      !allocate(real(8) :: rslt(size(buff,1), size(buff,2)))
  !end select
  !if (allocated(this%buffer)) then
    !rslt = this%buffer
  !else
    !call mpp_error(FATAL, 'get_buffer_2d: buffer not allocated')
  !endif
!end function get_buffer_2d

!> Gets real(r4_kind) buffer data from a 2d buffer object
!! called through this%get_buffer(r4_outputdata)
subroutine get_2d_real4 (this, buff_out)
  class(buffer2d_type), intent(in) :: this !< 2d allocated buffer object
  real(r4_kind), allocatable, intent(out)  :: buff_out(:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2)))
  select type (buff=>this%buffer)
    type is (real(r4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets real(r8_kind) buffer data from a 2d buffer object
!! called through this%get_buffer(r8_outputdata)
subroutine get_2d_real8 (this, buff_out)
  class(buffer2d_type), intent(in) :: this !< 2d allocated buffer object
  real(r8_kind), allocatable, intent(out)  :: buff_out(:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2)))
  select type (buff=>this%buffer)
    type is (real(r8_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i4_kind) buffer data from 2d buffer object
!> called through this%get_buffer(i4_outputdata)
subroutine get_2d_int4 (this, buff_out)
  class(buffer2d_type), intent(in) :: this !< 2d allocated buffer object
  integer(i4_kind), allocatable, intent(out)  :: buff_out(:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2)))
  select type (buff=>this%buffer)
    type is (integer(i4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i8_kind) buffer data from 2d buffer object
!> called through this%get_buffer(i8_outputdata)
subroutine get_2d_int8 (this, buff_out)
  class(buffer2d_type), intent(in) :: this !< 2d allocated buffer object
  integer(i8_kind), allocatable, intent(out)  :: buff_out(:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2)))
  select type (buff=>this%buffer)
    type is (integer(i8_kind))
      buff_out = buff
  end select
end subroutine


!> @brief Gets buffer data from buffer3d_type type
! @return copy of the buffer data
!function get_buffer_3d (this) &
!result(rslt)
  !class (buffer3d_type), intent(in) :: this !< 3D buffer object
  !class(*), allocatable :: rslt(:,:,:)
  !select type(buff=>this%buffer)
    !type is(integer(8))
      !allocate(integer(8) :: rslt(size(buff,1), size(buff,2), size(buff,3)))
    !type is(integer(4))
      !allocate(integer(4) :: rslt(size(buff,1), size(buff,2), size(buff,3)))
    !type is(real(4))
      !allocate(real(4) :: rslt(size(buff,1), size(buff,2), size(buff,3)))
    !type is(real(8))
      !allocate(real(8) :: rslt(size(buff,1), size(buff,2), size(buff,3)))
  !end select
  !if (allocated(this%buffer)) then
    !rslt = this%buffer
  !else
    !call mpp_error(FATAL, 'get_buffer_3d: buffer not allocated')
  !endif
!end function get_buffer_3d

!> Gets real(r4_kind) buffer data from a 3d buffer object
!! called through this%get_buffer(r4_outputdata)
subroutine get_3d_real4 (this, buff_out)
  class(buffer3d_type), intent(in) :: this !< 3d allocated buffer object
  real(r4_kind), allocatable, intent(out)  :: buff_out(:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3)))
  select type (buff=>this%buffer)
    type is (real(r4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets real(r8_kind) buffer data from a 3d buffer object
!! called through this%get_buffer(r8_outputdata)
subroutine get_3d_real8 (this, buff_out)
  class(buffer3d_type), intent(in) :: this !< 3d allocated buffer object
  real(r8_kind), allocatable, intent(out)  :: buff_out(:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3)))
  select type (buff=>this%buffer)
    type is (real(r8_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i4_kind) buffer data from 3d buffer object
!> called through this%get_buffer(i4_outputdata)
subroutine get_3d_int4 (this, buff_out)
  class(buffer3d_type), intent(in) :: this !< 3d allocated buffer object
  integer(i4_kind), allocatable, intent(out)  :: buff_out(:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3)))
  select type (buff=>this%buffer)
    type is (integer(i4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i8_kind) buffer data from 3d buffer object
!> called through this%get_buffer(i8_outputdata)
subroutine get_3d_int8 (this, buff_out)
  class(buffer3d_type), intent(in) :: this !< 3d allocated buffer object
  integer(i8_kind), allocatable, intent(out)  :: buff_out(:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3)))
  select type (buff=>this%buffer)
    type is (integer(i8_kind))
      buff_out = buff
  end select
end subroutine


!> @brief Gets buffer data from buffer4d_type type
!! @return copy of the buffer data
!function get_buffer_4d (this) &
!result(rslt)
  !class (buffer4d_type), intent(in) :: this !< 4D buffer object
  !class(*), allocatable :: rslt(:,:,:,:)
  !select type(buff=>this%buffer)
    !type is(integer(8))
      !allocate(integer(8) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4)))
    !type is(integer(4))
      !allocate(integer(4) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4)))
    !type is(real(4))
      !allocate(real(4) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4)))
    !type is(real(8))
      !allocate(real(8) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4)))
  !end select
  !if (allocated(this%buffer)) then
    !rslt = this%buffer
  !else
    !call mpp_error(FATAL, 'get_buffer_4d: buffer not allocated')
  !endif
!end function get_buffer_4d

!> Gets real(r4_kind) buffer data from a 4d buffer object
!! called through this%get_buffer(r4_outputdata)
subroutine get_4d_real4 (this, buff_out)
  class(buffer4d_type), intent(in) :: this !< 4d allocated buffer object
  real(r4_kind), allocatable, intent(out)  :: buff_out(:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), size(this%buffer,4)))
  select type (buff=>this%buffer)
    type is (real(r4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets real(r8_kind) buffer data from a 4d buffer object
!! called through this%get_buffer(r8_outputdata)
subroutine get_4d_real8 (this, buff_out)
  class(buffer4d_type), intent(in) :: this !< 4d allocated buffer object
  real(r8_kind), allocatable, intent(out)  :: buff_out(:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), size(this%buffer,4)))
  select type (buff=>this%buffer)
    type is (real(r8_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i4_kind) buffer data from 4d buffer object
!> called through this%get_buffer(i4_outputdata)
subroutine get_4d_int4 (this, buff_out)
  class(buffer4d_type), intent(in) :: this !< 4d allocated buffer object
  integer(i4_kind), allocatable, intent(out)  :: buff_out(:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), size(this%buffer,4)))
  select type (buff=>this%buffer)
    type is (integer(i4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i8_kind) buffer data from 4d buffer object
!> called through this%get_buffer(i8_outputdata)
subroutine get_4d_int8 (this, buff_out)
  class(buffer4d_type), intent(in) :: this !< 4d allocated buffer object
  integer(i8_kind), allocatable, intent(out)  :: buff_out(:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), size(this%buffer,4)))
  select type (buff=>this%buffer)
    type is (integer(i8_kind))
      buff_out = buff
  end select
end subroutine

!> @brief Gets buffer data from buffer5d_type type
!! @return copy of the buffer data
!function get_buffer_5d (this) &
!result(rslt)
  !class (buffer5d_type), intent(in) :: this !< 5D buffer object
  !class(*), allocatable :: rslt(:,:,:,:,:)
  !select type(buff=>this%buffer)
    !type is(integer(8))
      !allocate(integer(8) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4), size(buff,5)))
    !type is(integer(4))
      !allocate(integer(4) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4), size(buff,5)))
    !type is(real(4))
      !allocate(real(4) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4), size(buff,5)))
    !type is(real(8))
      !allocate(real(8) :: rslt(size(buff,1), size(buff,2), size(buff,3), size(buff,4), size(buff,5)))
  !end select
  !if (allocated(this%buffer)) then
    !rslt = this%buffer
  !else
    !call mpp_error(FATAL, 'get_buffer_5d: buffer not allocated')
  !endif
!end function get_buffer_5d

!> Gets real(r4_kind) buffer data from a 5d buffer object
!! called through this%get_buffer(r4_outputdata)
subroutine get_5d_real4 (this, buff_out)
  class(buffer5d_type), intent(in) :: this !< 5d allocated buffer object
  real(r4_kind), allocatable, intent(out)  :: buff_out(:,:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), &
                  & size(this%buffer,4), size(this%buffer,5)))
  select type (buff=>this%buffer)
    type is (real(r4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets real(r8_kind) buffer data from a 5d buffer object
!! called through this%get_buffer(r8_outputdata)
subroutine get_5d_real8 (this, buff_out)
  class(buffer5d_type), intent(in) :: this !< 5d allocated buffer object
  real(r8_kind), allocatable, intent(out)  :: buff_out(:,:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), &
                  & size(this%buffer,4), size(this%buffer,5)))
  select type (buff=>this%buffer)
    type is (real(r8_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i4_kind) buffer data from 5d buffer object
!> called through this%get_buffer(i4_outputdata)
subroutine get_5d_int4 (this, buff_out)
  class(buffer5d_type), intent(in) :: this !< 5d allocated buffer object
  integer(i4_kind), allocatable, intent(out)  :: buff_out(:,:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), &
                  & size(this%buffer,4), size(this%buffer,5)))
  select type (buff=>this%buffer)
    type is (integer(i4_kind))
      buff_out = buff
  end select
end subroutine

!> Gets integer(i8_kind) buffer data from 5d buffer object
!> called through this%get_buffer(i8_outputdata)
subroutine get_5d_int8 (this, buff_out)
  class(buffer5d_type), intent(in) :: this !< 5d allocated buffer object
  integer(i8_kind), allocatable, intent(out)  :: buff_out(:,:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  allocate(buff_out(size(this%buffer,1), size(this%buffer,2), size(this%buffer,3), &
                  & size(this%buffer,4), size(this%buffer)))
  select type (buff=>this%buffer)
    type is (integer(i8_kind))
      buff_out = buff
  end select
end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_0d (this, fillval)
  class(buffer0d_type), intent(inout) :: this !< scalar buffer object
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
  class(buffer1d_type), intent(inout) :: this !< 1D buffer object
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
  class(buffer2d_type), intent(inout) :: this !< 2D buffer object
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
  class(buffer3d_type), intent(inout) :: this !< 3D buffer object
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
  class(buffer4d_type), intent(inout) :: this !< allocated 4D buffer object
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
  class(buffer5d_type), intent(inout) :: this !< allocated 5D buffer object
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
  class(buffer0d_type), intent(inout) :: this !< allocated scalar buffer object
  class(*), intent(in)      :: input_data !< data to copy into buffer
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_1d: buffer not yet allocated')
  call this%initialize_buffer(input_data)
end subroutine add_to_buffer_0d

!> @brief Copy values ( from 1 to size(input_data)) into a 1d buffer object
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_1d(this, input_data)
  class(buffer1d_type), intent(inout) :: this !< allocated 1d buffer object
  class(*), intent(in)            :: input_data(:) !< data to copy into the buffer
  integer            :: n !< number of elements in input data
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  type_error = .false.
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
  class(buffer2d_type), intent(inout) :: this !< allocated 2d buffer object
  class(*), intent(in)             :: input_data(:,:) !< 2d data array to copy into buffer
  integer            :: n1, n2 !< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  type_error = .false.
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
  class(buffer3d_type), intent(inout) :: this !< allocated 3d buffer object
  class(*), intent(in)             :: input_data(:,:,:)!< 3d data array to copy into buffer
  integer            :: n1, n2, n3 !< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  type_error = .false.
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
  class(buffer4d_type), intent(inout) :: this !< allocated 4d buffer object
  class(*), intent(in)             :: input_data(:,:,:,:) !< 4d data to copy into buffer
  integer            :: n1, n2, n3, n4!< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  type_error = .false.
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
  class(buffer5d_type), intent(inout) :: this !< allocated 5d buffer object
  class(*), intent(in)             :: input_data(:,:,:,:,:) !< 5d data to copy into buffer
  integer            :: n1, n2, n3, n4, n5 !< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  type_error = .false.
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
