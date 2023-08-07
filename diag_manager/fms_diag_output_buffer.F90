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
module fms_diag_output_buffer_mod
#ifdef use_yaml
use platform_mod
use iso_c_binding
use time_manager_mod, only: time_type
use mpp_mod, only: mpp_error, FATAL
use diag_data_mod, only: DIAG_NULL, DIAG_NOT_REGISTERED, i4, i8, r4, r8
use fms2_io_mod, only: FmsNetcdfFile_t, write_data, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t
use fms_diag_yaml_mod, only: diag_yaml

implicit none

private

!> @brief Object that holds buffered data and other diagnostics
!! Abstract to ensure use through its extensions(buffer0-5d types)
type, abstract :: fmsDiagOutputBuffer_class
  integer, allocatable, private               :: buffer_id !< index in buffer list
  integer, allocatable, public :: num_elements(:) !< used in time-averaging
  class(*), allocatable, public :: count_0d(:) !< used in time-averaging along with
                                       !! counter which is stored in the child types (bufferNd)
  integer(i4_kind), public :: buffer_type !<set to allocated data type & kind value, one of i4,i8,r4,r8
  integer, allocatable, public :: buffer_dims(:) !< holds the size of each dimension in the buffer
  contains

  procedure :: flush_buffer
  procedure :: remap_buffer
  procedure :: set_buffer_id
  ! TODO could make these 'defered' ie. declared here but defined in each child type
  ! holding off cause the class(*) + polymorphism in here is probably already enough to upset the gods of compilation
  !procedure(allocate_buffer), deferred :: allocate_buffer
  !procedure, deferred :: get_buffer
  !procedure, deferred :: initialize_buffer

end type fmsDiagOutputBuffer_class

!> holds an allocated buffer0-5d object
type :: fmsDiagOutputBufferContainer_type
  class(fmsDiagOutputBuffer_class), allocatable :: diag_buffer_obj !< any 0-5d buffer object
  integer,                          allocatable :: axis_ids(:)     !< Axis ids for the buffer
  integer                                       :: field_id        !< The id of the field the buffer belongs to
  integer                                       :: yaml_id         !< The id of the yaml id the buffer belongs to

  contains
  procedure :: add_axis_ids
  procedure :: get_axis_ids
  procedure :: set_field_id
  procedure :: get_field_id
  procedure :: set_yaml_id
  procedure :: get_yaml_id
  procedure :: write_buffer
  !! These are needed because otherwise the write_data calls will go into the wrong interface
  procedure :: write_buffer_wrapper_netcdf
  procedure :: write_buffer_wrapper_domain
  procedure :: write_buffer_wrapper_u
end type

!> Scalar buffer type to extend fmsDiagBufferContainer_type
type, extends(fmsDiagOutputBuffer_class) :: outputBuffer0d_type
  class(*), allocatable :: buffer(:) !< "scalar" numeric buffer value
                                     !! will only be allocated to hold 1 value
  class(*), allocatable :: counter(:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_0d
  procedure :: initialize_buffer => initialize_buffer_0d
  procedure :: add_to_buffer => add_to_buffer_0d
  procedure :: get_buffer => get_0d

end type outputBuffer0d_type

!> 1D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagOutputBuffer_class) :: outputBuffer1d_type
  class(*), allocatable :: buffer(:) !< 1D numeric data array
  class(*), allocatable :: counter(:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_1d
  procedure :: initialize_buffer => initialize_buffer_1d
  procedure :: add_to_buffer => add_to_buffer_1d
  procedure :: get_buffer => get_1d
end type outputBuffer1d_type

!> 2D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagOutputBuffer_class) :: outputBuffer2d_type
  class(*), allocatable :: buffer(:,:) !< 2D numeric data array
  class(*), allocatable :: counter(:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_2d
  procedure :: initialize_buffer => initialize_buffer_2d
  procedure :: add_to_buffer => add_to_buffer_2d
  procedure :: get_buffer => get_2d
end type outputBuffer2d_type

!> 3D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagOutputBuffer_class) :: outputBuffer3d_type
  class(*), allocatable :: buffer(:,:,:) !< 3D numeric data array
  class(*), allocatable :: counter(:,:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_3d
  procedure :: initialize_buffer => initialize_buffer_3d
  procedure :: add_to_buffer => add_to_buffer_3d
  procedure :: get_buffer => get_3d
end type outputBuffer3d_type

!> 4D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagOutputBuffer_class) :: outputBuffer4d_type
  class(*), allocatable :: buffer(:,:,:,:) !< 4D numeric data array
  class(*), allocatable :: counter(:,:,:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_4d
  procedure :: initialize_buffer => initialize_buffer_4d
  procedure :: add_to_buffer => add_to_buffer_4d
  procedure :: get_buffer => get_4d
end type outputBuffer4d_type

!> 5D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagOutputBuffer_class) :: outputBuffer5d_type
  class(*), allocatable :: buffer(:,:,:,:,:) !< 5D numeric data array
  class(*), allocatable :: counter(:,:,:,:,:) !< (x,y,z, time-of-day) used in the time averaging functions
  contains
  procedure :: allocate_buffer => allocate_buffer_5d
  procedure :: initialize_buffer => initialize_buffer_5d
  procedure :: add_to_buffer => add_to_buffer_5d
  procedure :: get_buffer => get_5d
end type outputBuffer5d_type

! public types
public :: outputBuffer0d_type
public :: outputBuffer1d_type
public :: outputBuffer2d_type
public :: outputBuffer3d_type
public :: outputBuffer4d_type
public :: outputBuffer5d_type
public :: fmsDiagOutputBuffer_class
public :: fmsDiagOutputBufferContainer_type

! public routines
public :: fms_diag_output_buffer_init
public :: fms_diag_output_buffer_create_container

contains

!!--------module routines

!> Initializes a list of diag buffers
!> @returns true if allocation is successfull
logical function fms_diag_output_buffer_init(buffobjs, buff_list_size)
  type(fmsDiagOutputBufferContainer_type), allocatable, intent(out) :: buffobjs(:) !< an array of buffer container types
                                                                             !! to allocate
  integer, intent(in)                                          :: buff_list_size !< size of buffer array to allocate
  if (allocated(buffobjs)) call mpp_error(FATAL,'fms_diag_buffer_init: passed in buffobjs array is already allocated')
  allocate(buffobjs(buff_list_size))
  fms_diag_output_buffer_init = allocated(buffobjs)
end function fms_diag_output_buffer_init

!> Creates a container type encapsulating a new buffer object for the given dimensions.
!! The buffer object will still need to be allocated to a type via allocate_buffer() before use.
!> @result A fmsDiagBufferContainer_type that holds a bufferNd_type, where N is buff_dims
subroutine fms_diag_output_buffer_create_container(buff_dims, buffer_obj)
  integer,                                 intent(in)     :: buff_dims !< dimensions
  type(fmsDiagOutputBufferContainer_type), intent(inout)  :: buffer_obj

  character(len=5) :: dim_output !< string to output buff_dims on error

  select case (buff_dims)
    case (0)
      allocate(outputBuffer0d_type :: buffer_obj%diag_buffer_obj)
    case (1)
      allocate(outputBuffer1d_type :: buffer_obj%diag_buffer_obj)
    case (2)
      allocate(outputBuffer2d_type :: buffer_obj%diag_buffer_obj)
    case (3)
      allocate(outputBuffer3d_type :: buffer_obj%diag_buffer_obj)
    case (4)
      allocate(outputBuffer4d_type :: buffer_obj%diag_buffer_obj)
    case (5)
      allocate(outputBuffer5d_type :: buffer_obj%diag_buffer_obj)
    case default
      write( dim_output, *) buff_dims
      dim_output = adjustl(dim_output)
      call mpp_error(FATAL, 'fms_diag_buffer_create_container: invalid number of dimensions given:' // dim_output //&
                            '. Must be 0-5')
  end select
end subroutine fms_diag_output_buffer_create_container

!!--------generic routines for any fmsDiagBuffer_class objects

!> Setter for buffer_id for any buffer objects
subroutine set_buffer_id(this, id)
  class(fmsDiagOutputBuffer_class), intent(inout) :: this !< buffer object to set id for
  integer, intent(in)                       :: id !< positive integer id to set
  if (.not.allocated(this%buffer_id) ) allocate(this%buffer_id)
  this%buffer_id = id
end subroutine set_buffer_id

!> Remaps 0-5d data buffer from the given object onto a 5d array pointer.
!> @returns a 5D remapped buffer, with 1:1 for any added dimensions.
function remap_buffer(buffobj, field_name, has_diurnal_axis)
  class(fmsDiagOutputBuffer_class), target, intent(inout) :: buffobj !< any dimension buffer object
  class(*), pointer                                 :: remap_buffer(:,:,:,:,:)
  character(len=*), intent(in)                      :: field_name !< name of field for error output
  logical, intent(in) :: has_diurnal_axis !< true if the buffer has diurnal axis

  ! get num dimensions from type extension
  select type (buffobj)
    type is (outputBuffer0d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated" // &
                                                                 "for field:" // field_name)
      remap_buffer(1:size(buffobj%buffer,1), 1:1, 1:1, 1:1, 1:1) => buffobj%buffer
    type is (outputBuffer1d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated" // &
                                                                 "for field:" // field_name)
      remap_buffer(1:size(buffobj%buffer,1), 1:1, 1:1, 1:1, 1:1) => buffobj%buffer(1:size(buffobj%buffer,1))
    type is (outputBuffer2d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated" // &
                                                                 "for field:" // field_name)
      if (has_diurnal_axis) then
        remap_buffer(1:size(buffobj%buffer,1), 1:1, 1:1, 1:1, 1:size(buffobj%buffer,2)) => buffobj%buffer(:,:)
      else
        remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:1, 1:1, 1:1) => buffobj%buffer(:,:)
      end if
    type is (outputBuffer3d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated" // &
                                                                 "for field:" // field_name)
      if (has_diurnal_axis) then
        remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:1, 1:1, &
          1:size(buffobj%buffer,3)) => buffobj%buffer(:,:,:)
      else
        remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), &
          1:size(buffobj%buffer,3), 1:1, 1:1) => buffobj%buffer(:,:,:)
      end if
    type is (outputBuffer4d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated" // &
                                                                 "for field:" // field_name)
      if (has_diurnal_axis) then
        remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), &
                  1:1, 1:size(buffobj%buffer,4)) => buffobj%buffer(:,:,:,:)
      else
        remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), &
                  1:size(buffobj%buffer,4), 1:1) => buffobj%buffer(:,:,:,:)
      end if
    type is (outputBuffer5d_type)
      if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated" // &
                                                                 "for field:" // field_name)
      remap_buffer(1:size(buffobj%buffer,1), 1:size(buffobj%buffer,2), 1:size(buffobj%buffer,3), &
                   1:size(buffobj%buffer,4), 1:size(buffobj%buffer,5)) => buffobj%buffer(:,:,:,:,:)
    class default
      call mpp_error( FATAL, 'remap_buffer_pointer: invalid buffer type for remapping')
  end select

end function remap_buffer

!> Deallocates data fields from a buffer object.
subroutine flush_buffer(this)
  class(fmsDiagOutputBuffer_class), intent(inout) :: this !< any buffer object
  select type (this)
    type is (outputBuffer0d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
      if (allocated(this%counter)) deallocate(this%counter)
    type is (outputBuffer1d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
      if (allocated(this%counter)) deallocate(this%counter)
    type is (outputBuffer2d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
      if (allocated(this%counter)) deallocate(this%counter)
    type is (outputBuffer3d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
      if (allocated(this%counter)) deallocate(this%counter)
    type is (outputBuffer4d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
      if (allocated(this%counter)) deallocate(this%counter)
    type is (outputBuffer5d_type)
      if (allocated(this%buffer)) deallocate(this%buffer)
      if (allocated(this%counter)) deallocate(this%counter)
  end select
  if (allocated(this%buffer_id)) deallocate(this%buffer_id)
  if (allocated(this%count_0d)) deallocate(this%count_0d)
  if (allocated(this%num_elements)) deallocate(this%num_elements)
  if (allocated(this%buffer_dims)) deallocate(this%buffer_dims)
end subroutine flush_buffer

!! -----------Type-specific routines for buffer0-5d

!> Allocates scalar buffer data to the given buff_type.
subroutine allocate_buffer_0d(this, buff_type, field_name, diurnal_samples)
  class(outputBuffer0d_type), intent(inout), target :: this !< scalar buffer object
  class(*),intent(in) :: buff_type !< allocates to the given type, value does not matter
  character(len=*), intent(in) :: field_name !< field name for error output
  integer, intent(in),optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_0d: buffer already allocated for field:"// &
                                                    field_name)
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(1))
      allocate(integer(kind=i4_kind) :: this%counter(1))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0_i4_kind
      this%count_0d = 0_i4_kind
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(1))
      allocate(integer(kind=i8_kind) :: this%counter(1))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0_i8_kind
      this%count_0d = 0_i8_kind
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(1))
      allocate(real(kind=r4_kind) :: this%counter(1))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r4_kind
      this%count_0d = 0.0_r4_kind
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(1))
      allocate(real(kind=r8_kind) :: this%counter(1))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r8_kind
      this%count_0d = 0.0_r8_kind
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer_0d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4" // &
           "for field:" // field_name, &
           FATAL)
  end select

  allocate(this%num_elements(n_samples))
  allocate(this%buffer_dims(1))
  this%num_elements = 0
  this%buffer_dims(1) = 1

end subroutine allocate_buffer_0d

!> Allocates 1D buffer data to given buff_type.
subroutine allocate_buffer_1d(this, buff_type, buff_size, field_name, diurnal_samples)
  class(outputBuffer1d_type), intent(inout), target :: this !< scalar buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_size !< dimension bounds
  character(len=*), intent(in) :: field_name !< field name for error output
  integer, intent(in), optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_1d: buffer already allocated for field:" // &
                                                   field_name)
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_size))
      allocate(integer(kind=i4_kind) :: this%counter(buff_size))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0_i4_kind
      this%count_0d = 0_i4_kind
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_size))
      allocate(integer(kind=i8_kind) :: this%counter(buff_size))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0_i8_kind
      this%count_0d = 0_i8_kind
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_size))
      allocate(real(kind=r4_kind) :: this%count_0d(buff_size))
      allocate(real(kind=r4_kind) :: this%counter(n_samples))
      this%counter = 0.0_r4_kind
      this%count_0d = 0.0_r4_kind
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_size))
      allocate(real(kind=r8_kind) :: this%count_0d(buff_size))
      allocate(real(kind=r8_kind) :: this%counter(n_samples))
      this%counter = 0.0_r8_kind
      this%count_0d = 0.0_r8_kind
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer_1d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4 " // &
           "for field:" // field_name, &
           FATAL)
  end select

  allocate(this%num_elements(n_samples))
  allocate(this%buffer_dims(1))
  this%num_elements = 0
  this%count_0d   = 0
  this%buffer_dims(1) = buff_size

end subroutine allocate_buffer_1d

!> Allocates a 2D buffer to given buff_type.
subroutine allocate_buffer_2d(this, buff_type, buff_sizes, field_name, diurnal_samples)
  class(outputBuffer2d_type), intent(inout), target :: this !< 2D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(2) !< dimension sizes
  integer, intent(in),optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1
  character(len=*), intent(in) :: field_name !< field name for error output

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_2d: buffer already allocated for field: " // &
                                                    field_name)
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0_i4_kind
      this%count_0d = 0_i4_kind
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0_i8_kind
      this%count_0d = 0_i8_kind
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r4_kind
      this%count_0d = 0.0_r4_kind
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r8_kind) :: this%counter(buff_sizes(1), buff_sizes(2)))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r8_kind
      this%count_0d = 0.0_r8_kind
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer_1d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4" // &
           "for field:" // field_name, &
           FATAL)
  end select
  allocate(this%num_elements(n_samples))
  allocate(this%buffer_dims(2))
  this%num_elements = 0
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)

end subroutine allocate_buffer_2d

!> Allocates a 3D buffer to given buff_type.
subroutine allocate_buffer_3d(this, buff_type, buff_sizes, field_name, diurnal_samples)
  class(outputBuffer3d_type), intent(inout), target :: this !< 3D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(3) !< dimension sizes
  integer, intent(in),optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1
  character(len=*), intent(in) :: field_name !< field name for error output

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_3d: buffer already allocated for field" // &
                                                   field_name)
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0_i4_kind
      this%count_0d = 0_i4_kind
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0_i8_kind
      this%count_0d = 0_i8_kind
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r4_kind
      this%count_0d = 0.0_r4_kind
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r8_kind) :: this%counter( buff_sizes(1),buff_sizes(2), buff_sizes(3)))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r8_kind
      this%count_0d = 0.0_r8_kind
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer_3d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4" // &
           "for field:" // field_name, FATAL)
  end select

  allocate(this%num_elements(n_samples))
  this%num_elements = 0
  this%count_0d   = 0
  allocate(this%buffer_dims(3))
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)
  this%buffer_dims(3) = buff_sizes(3)

end subroutine allocate_buffer_3d

!> Allocates a 4D buffer to given buff_type.
subroutine allocate_buffer_4d(this, buff_type, buff_sizes, field_name, diurnal_samples)
  class(outputBuffer4d_type), intent(inout), target :: this !< 4D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(4) !< dimension buff_sizes
  character(len=*), intent(in) :: field_name !< field name for error output
  integer, intent(in),optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_4d: buffer already allocated for field:" // &
                                                   field_name)

  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0_i4_kind
      this%count_0d = 0_i4_kind
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0_i8_kind
      this%count_0d = 0_i8_kind
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r4_kind
      this%count_0d = 0.0_r4_kind
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r8_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r8_kind
      this%count_0d = 0.0_r8_kind
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer_4d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4" // &
           "for field:" // field_name, FATAL)
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

!> Allocates a 5D buffer to given buff_type.
subroutine allocate_buffer_5d(this, buff_type, buff_sizes, field_name, diurnal_samples)
  class(outputBuffer5d_type), intent(inout), target :: this !< 5D buffer object
  class(*),intent(in) :: buff_type !< allocates to the type of buff_type
  integer, intent(in) :: buff_sizes(5) !< dimension buff_sizes
  character(len=*), intent(in) :: field_name !< field name for error output
  integer, intent(in),optional :: diurnal_samples !< number of diurnal samples, passed in from diag_yaml
  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer_5d: buffer already allocated for field:" // &
                                                   field_name)
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                                  & buff_sizes(5)))
      allocate(integer(kind=i4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                   & buff_sizes(5)))
      allocate(integer(kind=i4_kind) :: this%count_0d(n_samples))
      this%counter = 0_i4_kind
      this%count_0d = 0_i4_kind
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                  & buff_sizes(5)))
      allocate(integer(kind=i8_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                  & buff_sizes(5)))
      allocate(integer(kind=i8_kind) :: this%count_0d(n_samples))
      this%counter = 0_i8_kind
      this%count_0d = 0_i8_kind
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                               & buff_sizes(5)))
      allocate(real(kind=r4_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                & buff_sizes(5)))
      allocate(real(kind=r4_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r4_kind
      this%count_0d = 0.0_r4_kind
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                               & buff_sizes(5)))
      allocate(real(kind=r8_kind) :: this%counter(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                & buff_sizes(5)))
      allocate(real(kind=r8_kind) :: this%count_0d(n_samples))
      this%counter = 0.0_r8_kind
      this%count_0d = 0.0_r8_kind
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer_5d", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4" // &
           "for field:" // field_name, FATAL)
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

!> Get routine for scalar buffers.
!! Sets the buff_out argument to the integer or real value currently stored in the buffer.
subroutine get_0d (this, buff_out, field_name)
  class(outputBuffer0d_type), intent(in) :: this !< 0d allocated buffer object
  class(*), allocatable, intent(out)  :: buff_out !< output of copied buffer data
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'get_0d(get_buffer): buffer not yet allocated for field:' &
                                                        & // field_name)
  select type (buff=>this%buffer)
    type is (real(r4_kind))
      allocate(real(r4_kind) :: buff_out)
      buff_out = buff(1)
    type is (real(r8_kind))
      allocate(real(r8_kind) :: buff_out)
      buff_out = buff(1)
    type is (integer(i4_kind))
      allocate(integer(i4_kind) :: buff_out)
      buff_out = buff(1)
    type is (integer(i8_kind))
      allocate(integer(i8_kind) :: buff_out)
      buff_out = buff(1)
    class default
      call mpp_error(FATAL, "get_0d: buffer allocated to invalid type(must be integer or real, kind size 4 or 8)." // &
                            field_name)
  end select
end subroutine

!> Get routine for 1D buffers.
!! Sets the buff_out argument to the integer or real array currently stored in the buffer.
subroutine get_1d (this, buff_out, field_name)
  class(outputBuffer1d_type), intent(in) :: this !< 1d allocated buffer object
  class(*), allocatable, intent(out)  :: buff_out(:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  integer(i4_kind) :: buff_size !< size for allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'get_0d(get_buffer): buffer not yet allocated for field:' &
                                                        & // field_name)
  buff_size = size(this%buffer,1)

  select type (buff=>this%buffer)
    type is (real(r4_kind))
      allocate(real(r4_kind) :: buff_out(buff_size))
      buff_out = buff
    type is (real(r8_kind))
      allocate(real(r8_kind) :: buff_out(buff_size))
      buff_out = buff
    type is (integer(i4_kind))
      allocate(integer(i4_kind) :: buff_out(buff_size))
      buff_out = buff
    type is (integer(i8_kind))
      allocate(integer(i8_kind) :: buff_out(buff_size))
      buff_out = buff
    class default
      call mpp_error(FATAL, "get_1d: buffer allocated to invalid type(must be integer or real, kind size 4 or 8)." // &
                            "field name: "// field_name)
  end select
end subroutine

!> Get routine for 2D buffers.
!! Sets the buff_out argument to the integer or real array currently stored in the buffer.
subroutine get_2d (this, buff_out, field_name)
  class(outputBuffer2d_type), intent(in) :: this !< 2d allocated buffer object
  class(*), allocatable, intent(out)  :: buff_out(:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  integer(i4_kind) :: buff_size(2) !< sizes for allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'get_2d(get_buffer): buffer not yet allocated for field:' &
                                                        & // field_name)
  buff_size(1) = size(this%buffer,1)
  buff_size(2) = size(this%buffer,2)

  select type (buff=>this%buffer)
    type is (real(r4_kind))
      allocate(real(r4_kind) :: buff_out(buff_size(1), buff_size(2)))
      buff_out = buff
    type is (real(r8_kind))
      allocate(real(r8_kind) :: buff_out(buff_size(1), buff_size(2)))
      buff_out = buff
    type is (integer(i4_kind))
      allocate(integer(i4_kind) :: buff_out(buff_size(1), buff_size(2)))
      buff_out = buff
    type is (integer(i8_kind))
      allocate(integer(i8_kind) :: buff_out(buff_size(1), buff_size(2)))
      buff_out = buff
    class default
      call mpp_error(FATAL, "get_2d: buffer allocated to invalid type(must be integer or real, kind size 4 or 8)." // &
                            "field name: "// field_name)

  end select
end subroutine

!> Get routine for 3D buffers.
!! Sets the buff_out argument to the integer or real array currently stored in the buffer.
subroutine get_3d (this, buff_out, field_name)
  class(outputBuffer3d_type), intent(in) :: this !< 3d allocated buffer object
  class(*), allocatable, intent(out)  :: buff_out(:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  integer(i4_kind) :: buff_size(3)!< sizes for allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'get_3d(get_buffer): buffer not yet allocated for field:' &
                                                        & // field_name)
  buff_size(1) = size(this%buffer,1)
  buff_size(2) = size(this%buffer,2)
  buff_size(3) = size(this%buffer,3)

  select type (buff=>this%buffer)
    type is (real(r4_kind))
      allocate(real(r4_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3)))
      buff_out = buff
    type is (real(r8_kind))
      allocate(real(r8_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3)))
      buff_out = buff
    type is (integer(i4_kind))
      allocate(integer(i4_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3)))
      buff_out = buff
    type is (integer(i8_kind))
      allocate(integer(i8_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3)))
      buff_out = buff
    class default
      call mpp_error(FATAL, "get_3d: buffer allocated to invalid type(must be integer or real, kind size 4 or 8)." // &
                            "field name: "// field_name)
  end select
end subroutine

!> Get routine for 4D buffers.
!! Sets the buff_out argument to the integer or real array currently stored in the buffer.
subroutine get_4d (this, buff_out, field_name)
  class(outputBuffer4d_type), intent(in) :: this !< 4d allocated buffer object
  class(*), allocatable, intent(out)  :: buff_out(:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  integer(i4_kind) :: buff_size(4)!< sizes for allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'get_4d(get_buffer): buffer not yet allocated for field:' &
                                                        & // field_name)
  buff_size(1) = size(this%buffer,1)
  buff_size(2) = size(this%buffer,2)
  buff_size(3) = size(this%buffer,3)
  buff_size(4) = size(this%buffer,4)

  select type (buff=>this%buffer)
    type is (real(r4_kind))
      allocate(real(r4_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4)))
      buff_out = buff
    type is (real(r8_kind))
      allocate(real(r8_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4)))
      buff_out = buff
    type is (integer(i4_kind))
      allocate(integer(i4_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4)))
      buff_out = buff
    type is (integer(i8_kind))
      allocate(integer(i8_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4)))
      buff_out = buff
    class default
      call mpp_error(FATAL, "get_4d: buffer allocated to invalid type(must be integer or real, kind size 4 or 8)." // &
                            "field name: "// field_name)
  end select
end subroutine

!> Get routine for 5D buffers.
!! Sets the buff_out argument to the integer or real array currently stored in the buffer.
subroutine get_5d (this, buff_out, field_name)
  class(outputBuffer5d_type), intent(in) :: this !< 5d allocated buffer object
  class(*), allocatable, intent(out)  :: buff_out(:,:,:,:,:) !< output of copied buffer data
                                                !! must be the same size as the allocated buffer
  integer(i4_kind) :: buff_size(5)!< sizes for allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'get_5d: buffer not yet allocated for field:' &
                                                        & // field_name)
  buff_size(1) = size(this%buffer,1)
  buff_size(2) = size(this%buffer,2)
  buff_size(3) = size(this%buffer,3)
  buff_size(4) = size(this%buffer,4)
  buff_size(5) = size(this%buffer,5)

  select type (buff=>this%buffer)
    type is (real(r4_kind))
      allocate(real(r4_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4), buff_size(5)))
      buff_out = buff
    type is (real(r8_kind))
      allocate(real(r8_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4), buff_size(5)))
      buff_out = buff
    type is (integer(i4_kind))
      allocate(integer(i4_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4), buff_size(5)))
      buff_out = buff
    type is (integer(i8_kind))
      allocate(integer(i8_kind) :: buff_out(buff_size(1), buff_size(2), buff_size(3), buff_size(4), buff_size(5)))
      buff_out = buff
    class default
      call mpp_error(FATAL, "get_5d: buffer allocated to invalid type(must be integer or real, kind size 4 or 8)." // &
                            "field name: "// field_name)
  end select
end subroutine

!> @brief Initializes a buffer to a given fill value.
subroutine initialize_buffer_0d (this, fillval, field_name)
  class(outputBuffer0d_type), intent(inout) :: this !< scalar buffer object
  class(*), intent(in) :: fillval !< fill value, must be same type as the allocated buffer in this
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_0d: field:'// field_name // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: fillval does not match up with allocated buffer type(r8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: fillval does not match up with allocated buffer type(r4_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: fillval does not match up with allocated buffer type(i8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_0d: fillval does not match up with allocated buffer type(i4_kind)' // &
                            ' for field' // field_name )
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_0d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine initialize_buffer_0d

!> @brief Initializes a buffer to a given fill value.
subroutine initialize_buffer_1d (this, fillval, field_name)
  class(outputBuffer1d_type), intent(inout) :: this !< 1D buffer object
  class(*), intent(in) :: fillval !< fill value, must be same type as the allocated buffer in this
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_1d: field:'// field_name // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: fillval does not match up with allocated buffer type(r8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: fillval does not match up with allocated buffer type(r4_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: fillval does not match up with allocated buffer type(i8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_1d: fillval does not match up with allocated buffer type(i4_kind)' // &
                            ' for field' // field_name )
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_1d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine initialize_buffer_1d

!> @brief Initializes a buffer to a given fill value.
subroutine initialize_buffer_2d (this, fillval, field_name)
  class(outputBuffer2d_type), intent(inout) :: this !< 2D buffer object
  class(*), intent(in) :: fillval !< fill value, must be same type as the allocated buffer in this
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_2d: field:'// field_name // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: fillval does not match up with allocated buffer type(r8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: fillval does not match up with allocated buffer type(r4_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: fillval does not match up with allocated buffer type(i8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_2d: fillval does not match up with allocated buffer type(i4_kind)' // &
                            ' for field' // field_name )
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_2d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine initialize_buffer_2d

!> @brief Initializes a buffer to a given fill value.
subroutine initialize_buffer_3d (this, fillval, field_name)
  class(outputBuffer3d_type), intent(inout) :: this !< 3D buffer object
  class(*), intent(in) :: fillval!< fill value, must be same type as the allocated buffer in this
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_3d: field:'// field_name // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: fillval does not match up with allocated buffer type(r8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: fillval does not match up with allocated buffer type(r4_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: fillval does not match up with allocated buffer type(i8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_3d: fillval does not match up with allocated buffer type(i4_kind)' // &
                            ' for field' // field_name )
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_3d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine initialize_buffer_3d

!> @brief Initializes a buffer to a given fill value.
subroutine initialize_buffer_4d (this, fillval, field_name)
  class(outputBuffer4d_type), intent(inout) :: this !< allocated 4D buffer object
  class(*), intent(in) :: fillval!< fill value, must be same type as the allocated buffer in this
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_4d: field:'// field_name // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: fillval does not match up with allocated buffer type(r8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: fillval does not match up with allocated buffer type(r4_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: fillval does not match up with allocated buffer type(i8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_4d: fillval does not match up with allocated buffer type(i4_kind)' // &
                            ' for field' // field_name )
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_4d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine initialize_buffer_4d

!> @brief Initializes a buffer to a given fill value.
subroutine initialize_buffer_5d (this, fillval, field_name)
  class(outputBuffer5d_type), intent(inout) :: this !< allocated 5D buffer object
  class(*), intent(in) :: fillval!< fill value, must be same type as the allocated buffer in this
  character(len=*), intent(in) :: field_name !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer_5d: field:'// field_name // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')
  ! have to check fill value and buffer types match
  select type(buff => this%buffer)
  type is(real(r8_kind))
    select type(fillval)
    type is(real(r8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: fillval does not match up with allocated buffer type(r8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(real(r4_kind))
    select type(fillval)
    type is(real(r4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: fillval does not match up with allocated buffer type(r4_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i8_kind))
    select type(fillval)
    type is(integer(i8_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: fillval does not match up with allocated buffer type(i8_kind)' // &
                            ' for field' // field_name )
    end select
  type is(integer(i4_kind))
    select type(fillval)
    type is(integer(i4_kind))
      buff = fillval
    class default
      call mpp_error(FATAL, 'initialize_buffer_5d: fillval does not match up with allocated buffer type(i4_kind)' // &
                            ' for field' // field_name )
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_5d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine initialize_buffer_5d

!> @brief Add values to 0d buffer.
!! This will just call the init routine since there's only one value.
!! @note input_data must match allocated type of buffer object.
subroutine add_to_buffer_0d(this, input_data, field_name)
  class(outputBuffer0d_type), intent(inout) :: this !< allocated scalar buffer object
  class(*), intent(in)      :: input_data !< data to copy into buffer
  character(len=*), intent(in) :: field_name !< field name for error output
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_0d: buffer not yet allocated for field:'// &
                                                           field_name)
  call this%initialize_buffer(input_data, field_name)
end subroutine add_to_buffer_0d

!> @brief Copy values ( from 1 to size(input_data)) into a 1d buffer object.
!! @note input_data must match allocated type of buffer object.
subroutine add_to_buffer_1d(this, input_data, field_name)
  class(outputBuffer1d_type), intent(inout) :: this !< allocated 1d buffer object
  class(*), intent(in)            :: input_data(:) !< data to copy into the buffer
  integer            :: n !< number of elements in input data
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output
  type_error = .false.
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_1d: buffer not yet allocated for field:' // &
                                                            field_name)
  n = SIZE(input_data)
  if( n .gt. SIZE(this%buffer)) call mpp_error( FATAL,"add_to_buffer_1d: input data larger than allocated buffer " // &
                                                      "for field: "// field_name)
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
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types'// &
                                         ' for field:' // field_name)
end subroutine add_to_buffer_1d

!> @brief Copy values ( from 1 to size(input_data)) into a 2d buffer object.
!! @note input_data must match allocated type of buffer object.
subroutine add_to_buffer_2d(this, input_data, field_name)
  class(outputBuffer2d_type), intent(inout) :: this !< allocated 2d buffer object
  class(*), intent(in)             :: input_data(:,:) !< 2d data array to copy into buffer
  integer            :: n1, n2 !< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output
  type_error = .false.
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_2d: buffer not yet allocated for field:' // &
                                                            field_name)
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2)) then
    call mpp_error( FATAL,"add_to_buffer_2d: input data larger than allocated buffer")
  endif
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
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types'//&
                                         ' for field:'// field_name)
end subroutine add_to_buffer_2d

!> @brief Copy values ( from 1 to size(input_data)) into a 3d buffer object.
!! @note input_data must match allocated type of buffer object.
subroutine add_to_buffer_3d(this, input_data, field_name)
  class(outputBuffer3d_type), intent(inout) :: this !< allocated 3d buffer object
  class(*), intent(in)             :: input_data(:,:,:)!< 3d data array to copy into buffer
  integer            :: n1, n2, n3 !< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output
  type_error = .false.
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_3d: buffer not yet allocated for field:'//&
                                                           field_name)
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  n3 = SIZE(input_data, 3)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2) .or. &
    n3 .gt. SIZE(this%buffer, 3)) then
    call mpp_error( FATAL,"add_to_buffer_3d: input data larger than allocated buffer for field:"//field_name)
  endif
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
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_3d: mismatch between allocated buffer and input data types'//&
                                         ' for field:'//field_name)
end subroutine add_to_buffer_3d

!> @brief Copy values ( from 1 to size(input_data)) into a 4d buffer object.
!! @note input_data must match allocated type of buffer object.
subroutine add_to_buffer_4d(this, input_data, field_name)
  class(outputBuffer4d_type), intent(inout) :: this !< allocated 4d buffer object
  class(*), intent(in)             :: input_data(:,:,:,:) !< 4d data to copy into buffer
  integer            :: n1, n2, n3, n4!< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output
  type_error = .false.
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_4d: buffer not yet allocated for field:'// &
                                                           field_name)
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  n3 = SIZE(input_data, 3)
  n4 = SIZE(input_data, 4)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2) .or. &
    n3 .gt. SIZE(this%buffer, 3) .or. n4 .gt. SIZE(this%buffer, 4)) then
    call mpp_error( FATAL,"add_to_buffer_4d: input data larger than allocated buffer for field:"//field_name)
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
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_4d: mismatch between allocated buffer and input data types'// &
                                         ' for field:' //field_name)
end subroutine add_to_buffer_4d

!> @brief Copy values (from 1 to size(input_data)) into a 5d buffer object.
!! @note input_data must match allocated type of buffer object.
subroutine add_to_buffer_5d(this, input_data, field_name)
  class(outputBuffer5d_type), intent(inout) :: this !< allocated 5d buffer object
  class(*), intent(in)             :: input_data(:,:,:,:,:) !< 5d data to copy into buffer
  integer            :: n1, n2, n3, n4, n5 !< number of elements per dimension
  logical            :: type_error !< set to true if mismatch between input_data and allocated buffer
  character(len=*), intent(in) :: field_name !< field name for error output
  type_error = .false.
  if( .not. allocated(this%buffer)) call mpp_error (FATAL, 'add_to_buffer_5d: buffer not yet allocated for field:'// &
                                                           field_name)
  n1 = SIZE(input_data, 1)
  n2 = SIZE(input_data, 2)
  n3 = SIZE(input_data, 3)
  n4 = SIZE(input_data, 4)
  n5 = SIZE(input_data, 5)
  if( n1 .gt. SIZE(this%buffer, 1) .or. n2 .gt. SIZE(this%buffer, 2) .or. &
    n3 .gt. SIZE(this%buffer, 3) .or. n4 .gt. SIZE(this%buffer, 4) .or. &
    n5 .gt. SIZE(this%buffer, 5)) then
    call mpp_error( FATAL,"add_to_buffer_4d: input data larger than allocated buffer for field:"//field_name)
  endif
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
  if( type_error ) call mpp_error (FATAL,'add_to_buffer_5d: mismatch between allocated buffer and input data types'//&
                                         'for field:'// field_name)
end subroutine add_to_buffer_5d

!> @brief Adds the axis ids to the buffer object
subroutine add_axis_ids(this, axis_ids)
  class(fmsDiagOutputBufferContainer_type), intent(inout) :: this        !< Buffer object
  integer,                                  intent(in)    :: axis_ids(:) !< Axis ids to add

  this%axis_ids = axis_ids
end subroutine

!> @brief Get the axis_ids for the buffer
!! @return Axis_ids, if the buffer doesn't have axis ids it returns diag_null
function get_axis_ids(this) &
  result(res)

  class(fmsDiagOutputBufferContainer_type), intent(inout) :: this        !< Buffer object
  integer, allocatable :: res(:)

  if (allocated(this%axis_ids)) then
    res = this%axis_ids
  else
    allocate(res(1))
    res = diag_null
  endif
end function

!> @brief Get the field id of the buffer
!! @return the field id of the buffer
function get_field_id(this) &
  result(res)

  class(fmsDiagOutputBufferContainer_type), intent(in) :: this        !< Buffer object
  integer :: res

  res = this%field_id
end function get_field_id

!> @brief set the field id of the buffer
subroutine set_field_id(this, field_id)
  class(fmsDiagOutputBufferContainer_type), intent(inout) :: this        !< Buffer object
  integer,                                  intent(in)    :: field_id    !< field id of the buffer

  this%field_id = field_id
end subroutine set_field_id

!> @brief set the field id of the buffer
subroutine set_yaml_id(this, yaml_id)
  class(fmsDiagOutputBufferContainer_type), intent(inout) :: this        !< Buffer object
  integer,                                  intent(in)    :: yaml_id     !< yaml id of the buffer

  this%yaml_id = yaml_id
end subroutine set_yaml_id

!> @brief Get the yaml id of the buffer
!! @return the yaml id of the buffer
function get_yaml_id(this) &
  result(res)

  class(fmsDiagOutputBufferContainer_type), intent(in) :: this        !< Buffer object
  integer :: res

  res = this%yaml_id
end function get_yaml_id

!> @brief Write the buffer to the file
subroutine write_buffer(this, fileobj, unlim_dim_level)
  class(fmsDiagOutputBufferContainer_type), intent(in) :: this            !< buffer object to write
  class(FmsNetcdfFile_t),                   intent(in) :: fileobj         !< fileobj to write to
  integer, optional,                        intent(in) :: unlim_dim_level !< unlimited dimension

  select type(fileobj)
  type is (FmsNetcdfFile_t)
    call this%write_buffer_wrapper_netcdf(fileobj, unlim_dim_level=unlim_dim_level)
  type is (FmsNetcdfDomainFile_t)
    call this%write_buffer_wrapper_domain(fileobj, unlim_dim_level=unlim_dim_level)
  type is (FmsNetcdfUnstructuredDomainFile_t)
    call this%write_buffer_wrapper_u(fileobj, unlim_dim_level=unlim_dim_level)
  end select
end subroutine write_buffer

!> @brief Write the buffer to the FmsNetcdfFile_t fileobj
subroutine write_buffer_wrapper_netcdf(this, fileobj, unlim_dim_level)
  class(fmsDiagOutputBufferContainer_type), intent(in) :: this            !< buffer object to write
  type(FmsNetcdfFile_t),                    intent(in) :: fileobj         !< fileobj to write to
  integer, optional,                        intent(in) :: unlim_dim_level !< unlimited dimension

  character(len=:), allocatable :: varname !< name of the variable

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select type(buffer_obj=>this%diag_buffer_obj)
  type is (outputBuffer0d_type)
    call write_data(fileobj, varname, buffer_obj%buffer(1), unlim_dim_level=unlim_dim_level)
  type is (outputBuffer1d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer2d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer3d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer4d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer5d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  end select
end subroutine

!> @brief Write the buffer to the FmsNetcdfDomainFile_t fileobj
subroutine write_buffer_wrapper_domain(this, fileobj, unlim_dim_level)
  class(fmsDiagOutputBufferContainer_type), intent(in) :: this            !< buffer object to write
  type(FmsNetcdfDomainFile_t),              intent(in) :: fileobj         !< fileobj to write to
  integer, optional,                        intent(in) :: unlim_dim_level !< unlimited dimension

  character(len=:), allocatable :: varname !< name of the variable

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select type(buffer_obj=>this%diag_buffer_obj)
  type is (outputBuffer0d_type)
    call write_data(fileobj, varname, buffer_obj%buffer(1), unlim_dim_level=unlim_dim_level)
  type is (outputBuffer1d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer2d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer3d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer4d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer5d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  end select
end subroutine

!> @brief Write the buffer to the FmsNetcdfUnstructuredDomainFile_t fileobj
subroutine write_buffer_wrapper_u(this, fileobj, unlim_dim_level)
  class(fmsDiagOutputBufferContainer_type), intent(in) :: this            !< buffer object to write
  type(FmsNetcdfUnstructuredDomainFile_t),  intent(in) :: fileobj         !< fileobj to write to
  integer, optional,                        intent(in) :: unlim_dim_level !< unlimited dimension

  character(len=:), allocatable :: varname !< name of the variable

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select type(buffer_obj=>this%diag_buffer_obj)
  type is (outputBuffer0d_type)
    call write_data(fileobj, varname, buffer_obj%buffer(1), unlim_dim_level=unlim_dim_level)
  type is (outputBuffer1d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer2d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer3d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer4d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  type is (outputBuffer5d_type)
    call write_data(fileobj, varname, buffer_obj%buffer, unlim_dim_level=unlim_dim_level)
  end select
end subroutine
#endif
end module fms_diag_output_buffer_mod
