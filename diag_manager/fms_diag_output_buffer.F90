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
use time_manager_mod, only: time_type, operator(==), operator(>)
use mpp_mod, only: mpp_error, FATAL, NOTE
use diag_data_mod, only: DIAG_NULL, DIAG_NOT_REGISTERED, i4, i8, r4, r8, get_base_time, MIN_VALUE, MAX_VALUE, EMPTY, &
                         time_min, time_max
use fms2_io_mod, only: FmsNetcdfFile_t, write_data, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t
use fms_diag_yaml_mod, only: diag_yaml
use fms_diag_bbox_mod, only: fmsDiagIbounds_type
use fms_diag_reduction_methods_mod, only: do_time_none, do_time_min, do_time_max, do_time_sum_update, time_update_done
use fms_diag_time_utils_mod, only: diag_time_inc

implicit none

private

!> holds an allocated buffer0-5d object
type :: fmsDiagOutputBuffer_type
  integer               :: buffer_id          !< index in buffer list
  integer(i4_kind)      :: buffer_type        !< set to allocated data type & kind value, one of i4,i8,r4,r8
  class(*), allocatable :: buffer(:,:,:,:,:)  !< 5D numeric data array
  integer               :: ndim               !< Number of dimensions for each variable
  integer,  allocatable :: buffer_dims(:)     !< holds the size of each dimension in the buffer
  real(r8_kind)         :: weight_sum !< (x,y,z, time-of-day) used in the time averaging functions
  integer,  allocatable :: num_elements(:)    !< used in time-averaging
  integer,  allocatable :: axis_ids(:)        !< Axis ids for the buffer
  integer               :: field_id           !< The id of the field the buffer belongs to
  integer               :: yaml_id            !< The id of the yaml id the buffer belongs to
  logical               :: done_with_math     !< .True. if done doing the math
  type(time_type)       :: time               !< The last time the data was received

  contains
  procedure :: add_axis_ids
  procedure :: get_axis_ids
  procedure :: set_field_id
  procedure :: get_field_id
  procedure :: set_yaml_id
  procedure :: get_yaml_id
  procedure :: init_buffer_time
  procedure :: update_buffer_time
  procedure :: is_done_with_math
  procedure :: set_done_with_math
  procedure :: write_buffer
  !! These are needed because otherwise the write_data calls will go into the wrong interface
  procedure :: write_buffer_wrapper_netcdf
  procedure :: write_buffer_wrapper_domain
  procedure :: write_buffer_wrapper_u
  procedure :: allocate_buffer
  procedure :: initialize_buffer
  procedure :: get_buffer
  procedure :: flush_buffer
  procedure :: do_time_none_wrapper
  procedure :: do_time_min_wrapper
  procedure :: do_time_max_wrapper
  procedure :: do_time_sum_wrapper
  procedure :: diag_reduction_done_wrapper
  procedure :: get_buffer_dims
end type fmsDiagOutputBuffer_type

! public types
public :: fmsDiagOutputBuffer_type

! public routines
public :: fms_diag_output_buffer_init

contains

!!--------module routines

!> Initializes a list of diag buffers
!> @returns true if allocation is successfull
logical function fms_diag_output_buffer_init(buffobjs, buff_list_size)
  type(fmsDiagOutputBuffer_type), allocatable, intent(out) :: buffobjs(:)    !< an array of buffer container types
                                                                        !! to allocate
  integer,                                     intent(in)  :: buff_list_size !< size of buffer array to allocate

  if (allocated(buffobjs)) call mpp_error(FATAL,'fms_diag_buffer_init: passed in buffobjs array is already allocated')
  allocate(buffobjs(buff_list_size))
  fms_diag_output_buffer_init = allocated(buffobjs)
end function fms_diag_output_buffer_init

!!--------generic routines for any fmsDiagBuffer_class objects

!> Setter for buffer_id for any buffer objects
subroutine set_buffer_id(this, id)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this !< buffer object to set id for
  integer,                         intent(in)    :: id   !< positive integer id to set

  this%buffer_id = id
end subroutine set_buffer_id

!> Deallocates data fields from a buffer object.
subroutine flush_buffer(this)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this !< any buffer object

  this%buffer_id   = diag_null
  this%buffer_type = diag_null
  this%ndim        = diag_null
  this%field_id    = diag_null
  this%yaml_id     = diag_null
  if (allocated(this%buffer))       deallocate(this%buffer)
  if (allocated(this%buffer_dims))  deallocate(this%buffer_dims)
  if (allocated(this%num_elements)) deallocate(this%num_elements)
  if (allocated(this%axis_ids))     deallocate(this%axis_ids)
end subroutine flush_buffer

!> Allocates a 5D buffer to given buff_type.
subroutine allocate_buffer(this, buff_type, ndim, buff_sizes, field_name, diurnal_samples)
  class(fmsDiagOutputBuffer_type), intent(inout), target :: this            !< 5D buffer object
  class(*),                        intent(in)            :: buff_type       !< allocates to the type of buff_type
  integer,                         intent(in)            :: ndim            !< Number of dimension
  integer,                         intent(in)            :: buff_sizes(5)   !< dimension buff_sizes
  character(len=*),                intent(in)            :: field_name      !< field name for error output
  integer, optional,               intent(in)            :: diurnal_samples !< number of diurnal samples

  integer :: n_samples !< number of diurnal samples, defaults to 1

  if(present(diurnal_samples)) then
    n_samples = diurnal_samples
  else
    n_samples = 1
  endif

  this%ndim =ndim
  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer: buffer already allocated for field:" // &
                                                   field_name)
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                                  & buff_sizes(5)))
      this%weight_sum = 0.0_r4_kind
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                  & buff_sizes(5)))
      this%weight_sum = 0.0_r8_kind
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                               & buff_sizes(5)))
      this%weight_sum = 0.0_r4_kind
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                               & buff_sizes(5)))
      this%weight_sum = 0.0_r8_kind
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4" // &
           "for field:" // field_name, FATAL)
  end select
  allocate(this%num_elements(n_samples))
  this%num_elements = 0
  this%done_with_math = .false.
  allocate(this%buffer_dims(5))
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)
  this%buffer_dims(3) = buff_sizes(3)
  this%buffer_dims(4) = buff_sizes(4)
  this%buffer_dims(5) = buff_sizes(5)
end subroutine allocate_buffer

!> Get routine for 5D buffers.
!! Sets the buff_out argument to the integer or real array currently stored in the buffer.
subroutine get_buffer (this, buff_out, field_name)
  class(fmsDiagOutputBuffer_type), intent(in)   :: this                !< 5d allocated buffer object
  class(*), allocatable,           intent(out)  :: buff_out(:,:,:,:,:) !< output of copied buffer data
                                                                       !! must be the same size as the allocated buffer
  character(len=*),                intent(in)   :: field_name          !< field name for error output

  integer(i4_kind) :: buff_size(5)!< sizes for allocated buffer

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'get_buffer: buffer not yet allocated for field:' &
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
      call mpp_error(FATAL, "get_buffer: buffer allocated to invalid type(must be integer or real, kind size 4 or 8)."&
                            //"field name: "// field_name)
  end select
end subroutine

!> @brief Initializes a buffer based on the reduction method
subroutine initialize_buffer (this, reduction_method, field_name)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this             !< allocated 5D buffer object
  integer,                         intent(in)    :: reduction_method !< The reduction method for the field
  character(len=*),                intent(in)    :: field_name       !< field name for error output

  if(.not. allocated(this%buffer)) call mpp_error(FATAL, 'initialize_buffer: field:'// field_name // &
      'buffer not yet allocated, allocate_buffer() must be called on this object first.')

  select type(buff => this%buffer)
  type is(real(r8_kind))
    select case (reduction_method)
    case (time_min)
      buff = real(MIN_VALUE, kind=r8_kind)
    case (time_max)
      buff = real(MAX_VALUE, kind=r8_kind)
    case default
      buff = real(EMPTY, kind=r8_kind)
    end select
  type is(real(r4_kind))
    select case (reduction_method)
    case (time_min)
      buff = real(MIN_VALUE, kind=r4_kind)
    case (time_max)
      buff = real(MAX_VALUE, kind=r4_kind)
    case default
      buff = real(EMPTY, kind=r4_kind)
    end select
  type is(integer(i8_kind))
    select case (reduction_method)
    case (time_min)
      buff = int(MIN_VALUE, kind=i8_kind)
    case (time_max)
      buff = int(MAX_VALUE, kind=i8_kind)
    case default
      buff = int(EMPTY, kind=i8_kind)
    end select
  type is(integer(i4_kind))
    select case (reduction_method)
    case (time_min)
      buff = int(MIN_VALUE, kind=i4_kind)
    case (time_max)
      buff = int(MAX_VALUE, kind=i4_kind)
    case default
      buff = int(EMPTY, kind=i4_kind)
    end select
  class default
    call mpp_error(FATAL, 'initialize buffer_5d: buffer allocated to invalid data type, this shouldnt happen')
  end select

end subroutine initialize_buffer

!> @brief Adds the axis ids to the buffer object
subroutine add_axis_ids(this, axis_ids)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  integer,                         intent(in)    :: axis_ids(:) !< Axis ids to add

  this%axis_ids = axis_ids
end subroutine

!> @brief Get the axis_ids for the buffer
!! @return Axis_ids, if the buffer doesn't have axis ids it returns diag_null
function get_axis_ids(this) &
  result(res)

  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
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
  class(fmsDiagOutputBuffer_type), intent(in) :: this        !< Buffer object
  integer :: res

  res = this%field_id
end function get_field_id

!> @brief set the field id of the buffer
subroutine set_field_id(this, field_id)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  integer,                         intent(in)    :: field_id    !< field id of the buffer

  this%field_id = field_id
end subroutine set_field_id

!> @brief set the field id of the buffer
subroutine set_yaml_id(this, yaml_id)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  integer,                         intent(in)    :: yaml_id     !< yaml id of the buffer

  this%yaml_id = yaml_id
end subroutine set_yaml_id

!> @brief inits the buffer time for the buffer
subroutine init_buffer_time(this, time)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  type(time_type), optional,       intent(in)    :: time        !< time to add to the buffer

  if (present(time)) then
    this%time = time
  else
    this%time = get_base_time()
  endif
end subroutine init_buffer_time

!> @brief Update the buffer time if it is a new time
!! @return .true. if the buffer was updated
function update_buffer_time(this, time) &
  result(res)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  type(time_type),                 intent(in)    :: time        !< time to add to the buffer

  logical :: res

  if (time > this%time) then
    this%time = time
    res = .true.
  else
    res = .false.
  endif
end function

!> @brief Determine if finished with math
!! @return this%done_with_math
function is_done_with_math(this) &
  result(res)
  class(fmsDiagOutputBuffer_type), intent(in) :: this        !< Buffer object
  logical :: res

  res = this%done_with_math
end function is_done_with_math

!> @brief Set done_with_math to .true.
subroutine set_done_with_math(this)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  integer :: res

  this%done_with_math = .true.
end subroutine set_done_with_math

!> @brief Get the yaml id of the buffer
!! @return the yaml id of the buffer
function get_yaml_id(this) &
  result(res)

  class(fmsDiagOutputBuffer_type), intent(in) :: this        !< Buffer object
  integer :: res

  res = this%yaml_id
end function get_yaml_id

!> @brief Write the buffer to the file
subroutine write_buffer(this, fms2io_fileobj, unlim_dim_level)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this            !< buffer object to write
  class(FmsNetcdfFile_t),          intent(in)    :: fms2io_fileobj  !< fileobj to write to
  integer, optional,               intent(in)    :: unlim_dim_level !< unlimited dimension

  select type(fms2io_fileobj)
  type is (FmsNetcdfFile_t)
    call this%write_buffer_wrapper_netcdf(fms2io_fileobj, unlim_dim_level=unlim_dim_level)
  type is (FmsNetcdfDomainFile_t)
    call this%write_buffer_wrapper_domain(fms2io_fileobj, unlim_dim_level=unlim_dim_level)
  type is (FmsNetcdfUnstructuredDomainFile_t)
    call this%write_buffer_wrapper_u(fms2io_fileobj, unlim_dim_level=unlim_dim_level)
  class default
    call mpp_error(FATAL, "The file "//trim(fms2io_fileobj%path)//" is not one of the accepted types"//&
      " only FmsNetcdfFile_t, FmsNetcdfDomainFile_t, and FmsNetcdfUnstructuredDomainFile_t are accepted.")
  end select

  call this%initialize_buffer(diag_yaml%diag_fields(this%yaml_id)%get_var_reduction(), &
    diag_yaml%diag_fields(this%yaml_id)%get_var_outname())
  !TODO Set the counters back to 0
end subroutine write_buffer

!> @brief Write the buffer to the FmsNetcdfFile_t fms2io_fileobj
subroutine write_buffer_wrapper_netcdf(this, fms2io_fileobj, unlim_dim_level)
  class(fmsDiagOutputBuffer_type),  intent(in) :: this            !< buffer object to write
  type(FmsNetcdfFile_t),            intent(in) :: fms2io_fileobj  !< fileobj to write to
  integer, optional,                intent(in) :: unlim_dim_level !< unlimited dimension

  character(len=:), allocatable :: varname !< name of the variable

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select case(this%ndim)
  case (0)
    call write_data(fms2io_fileobj, varname, this%buffer(1,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (1)
    call write_data(fms2io_fileobj, varname, this%buffer(:,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (2)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,1,1,1), unlim_dim_level=unlim_dim_level)
  case (3)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,1,1), unlim_dim_level=unlim_dim_level)
  case (4)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,:,1), unlim_dim_level=unlim_dim_level)
  case (5)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,:,:), unlim_dim_level=unlim_dim_level)
  end select
end subroutine write_buffer_wrapper_netcdf

!> @brief Write the buffer to the FmsNetcdfDomainFile_t fms2io_fileobj
subroutine write_buffer_wrapper_domain(this, fms2io_fileobj, unlim_dim_level)
  class(fmsDiagOutputBuffer_type),    intent(in) :: this            !< buffer object to write
  type(FmsNetcdfDomainFile_t),        intent(in) :: fms2io_fileobj  !< fileobj to write to
  integer, optional,                  intent(in) :: unlim_dim_level !< unlimited dimension

  character(len=:), allocatable :: varname !< name of the variable

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select case(this%ndim)
  case (0)
    call write_data(fms2io_fileobj, varname, this%buffer(1,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (1)
    call write_data(fms2io_fileobj, varname, this%buffer(:,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (2)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,1,1,1), unlim_dim_level=unlim_dim_level)
  case (3)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,1,1), unlim_dim_level=unlim_dim_level)
  case (4)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,:,1), unlim_dim_level=unlim_dim_level)
  case (5)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,:,:), unlim_dim_level=unlim_dim_level)
  end select
end subroutine write_buffer_wrapper_domain

!> @brief Write the buffer to the FmsNetcdfUnstructuredDomainFile_t fms2io_fileobj
subroutine write_buffer_wrapper_u(this, fms2io_fileobj, unlim_dim_level)
  class(fmsDiagOutputBuffer_type),                 intent(in) :: this            !< buffer object to write
  type(FmsNetcdfUnstructuredDomainFile_t),         intent(in) :: fms2io_fileobj  !< fileobj to write to
  integer, optional,                               intent(in) :: unlim_dim_level !< unlimited dimension

  character(len=:), allocatable :: varname !< name of the variable

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select case(this%ndim)
  case (0)
    call write_data(fms2io_fileobj, varname, this%buffer(1,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (1)
    call write_data(fms2io_fileobj, varname, this%buffer(:,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (2)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,1,1,1), unlim_dim_level=unlim_dim_level)
  case (3)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,1,1), unlim_dim_level=unlim_dim_level)
  case (4)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,:,1), unlim_dim_level=unlim_dim_level)
  case (5)
    call write_data(fms2io_fileobj, varname, this%buffer(:,:,:,:,:), unlim_dim_level=unlim_dim_level)
  end select
end subroutine write_buffer_wrapper_u

!> @brief Does the time_none reduction method on the buffer object
!! @return Error message if the math was not successful
function do_time_none_wrapper(this, field_data, mask, is_masked, bounds_in, bounds_out, missing_value) &
  result(err_msg)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this                !< buffer object to write
  class(*),                        intent(in)    :: field_data(:,:,:,:) !< Buffer data for current time
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_in           !< Indicies for the buffer passed in
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_out          !< Indicies for the output buffer
  logical,                         intent(in)    :: mask(:,:,:,:)       !< Mask for the field
  logical,                         intent(in)    :: is_masked           !< .True. if the field has a mask
  real(kind=r8_kind),              intent(in)    :: missing_value       !< Missing_value for data points that are masked
  character(len=50) :: err_msg

  !TODO This will be expanded for integers
  err_msg = ""
  select type (output_buffer => this%buffer)
    type is (real(kind=r8_kind))
      select type (field_data)
      type is (real(kind=r8_kind))
        call do_time_none(output_buffer, field_data, mask, is_masked, bounds_in, bounds_out, missing_value)
      class default
        err_msg="do_time_none_wrapper::the output buffer and the buffer send in are not of the same type (r8_kind)"
      end select
    type is (real(kind=r4_kind))
      select type (field_data)
      type is (real(kind=r4_kind))
        call do_time_none(output_buffer, field_data, mask, is_masked, bounds_in, bounds_out, &
          real(missing_value, kind=r4_kind))
      class default
        err_msg="do_time_none_wrapper::the output buffer and the buffer send in are not of the same type (r4_kind)"
      end select
  end select
end function do_time_none_wrapper

!> @brief Does the time_min reduction method on the buffer object
!! @return Error message if the math was not successful
function do_time_min_wrapper(this, field_data, mask, is_masked, bounds_in, bounds_out, missing_value) &
  result(err_msg)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this                !< buffer object to write
  class(*),                        intent(in)    :: field_data(:,:,:,:) !< Buffer data for current time
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_in           !< Indicies for the buffer passed in
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_out          !< Indicies for the output buffer
  logical,                         intent(in)    :: mask(:,:,:,:)       !< Mask for the field
  logical,                         intent(in)    :: is_masked           !< .True. if the field has a mask
  real(kind=r8_kind),              intent(in)    :: missing_value       !< Missing_value for data points that are masked
  character(len=50) :: err_msg

  !TODO This will be expanded for integers
  err_msg = ""
  select type (output_buffer => this%buffer)
    type is (real(kind=r8_kind))
      select type (field_data)
      type is (real(kind=r8_kind))
        call do_time_min(output_buffer, field_data, mask, is_masked, bounds_in, bounds_out, missing_value)
      class default
        err_msg="do_time_min_wrapper::the output buffer and the buffer send in are not of the same type (r8_kind)"
      end select
    type is (real(kind=r4_kind))
      select type (field_data)
      type is (real(kind=r4_kind))
        call do_time_min(output_buffer, field_data, mask, is_masked, bounds_in, bounds_out, &
          real(missing_value, kind=r4_kind))
      class default
        err_msg="do_time_min_wrapper::the output buffer and the buffer send in are not of the same type (r4_kind)"
      end select
  end select
end function do_time_min_wrapper

!> @brief Does the time_min reduction method on the buffer object
!! @return Error message if the math was not successful
function do_time_max_wrapper(this, field_data, mask, is_masked, bounds_in, bounds_out, missing_value) &
  result(err_msg)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this                !< buffer object to write
  class(*),                        intent(in)    :: field_data(:,:,:,:) !< Buffer data for current time
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_in           !< Indicies for the buffer passed in
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_out          !< Indicies for the output buffer
  logical,                         intent(in)    :: mask(:,:,:,:)       !< Mask for the field
  logical,                         intent(in)    :: is_masked           !< .True. if the field has a mask
  real(kind=r8_kind),              intent(in)    :: missing_value       !< Missing_value for data points that are masked
  character(len=50) :: err_msg

  !TODO This will be expanded for integers
  err_msg = ""
  select type (output_buffer => this%buffer)
    type is (real(kind=r8_kind))
      select type (field_data)
      type is (real(kind=r8_kind))
        call do_time_max(output_buffer, field_data, mask, is_masked, bounds_in, bounds_out, missing_value)
      class default
        err_msg="do_time_max_wrapper::the output buffer and the buffer send in are not of the same type (r8_kind)"
      end select
    type is (real(kind=r4_kind))
      select type (field_data)
      type is (real(kind=r4_kind))
        call do_time_max(output_buffer, field_data, mask, is_masked, bounds_in, bounds_out, &
          real(missing_value, kind=r4_kind))
      class default
        err_msg="do_time_max_wrapper::the output buffer and the buffer send in are not of the same type (r4_kind)"
      end select
  end select
end function do_time_max_wrapper

!> @brief Does the time_sum reduction method on the buffer object
!! @return Error message if the math was not successful
function do_time_sum_wrapper(this, field_data, mask, is_masked, bounds_in, bounds_out, missing_value, &
                             increase_counter) &
  result(err_msg)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this                !< buffer object to write
  class(*),                        intent(in)    :: field_data(:,:,:,:) !< Buffer data for current time
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_in           !< Indicies for the buffer passed in
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_out          !< Indicies for the output buffer
  logical,                         intent(in)    :: mask(:,:,:,:)       !< Mask for the field
  logical,                         intent(in)    :: is_masked           !< .True. if the field has a mask
  real(kind=r8_kind),              intent(in)    :: missing_value       !< Missing_value for data points that are masked
  logical,                         intent(in)    :: increase_counter    !< .True. if data has not been received for
                                                                        !! time, so the counter needs to be increased
  character(len=50) :: err_msg

  !TODO This will be expanded for integers
  err_msg = ""
  select type (output_buffer => this%buffer)
    type is (real(kind=r8_kind))
      select type (field_data)
      type is (real(kind=r8_kind))
        call do_time_sum_update(output_buffer, this%weight_sum, field_data, mask, is_masked, &
                                bounds_in, bounds_out, missing_value, increase_counter)
      class default
        err_msg="do_time_sum_wrapper::the output buffer and the buffer send in are not of the same type (r8_kind)"
      end select
    type is (real(kind=r4_kind))
      select type (field_data)
      type is (real(kind=r4_kind))
        call do_time_sum_update(output_buffer, this%weight_sum, field_data, mask, is_masked, bounds_in, bounds_out, &
          real(missing_value, kind=r4_kind), increase_counter)
      class default
        err_msg="do_time_sum_wrapper::the output buffer and the buffer send in are not of the same type (r4_kind)"
      end select
    class default
      err_msg="do_time_sum_wrapper::the output buffer is not a valid type, must be real(r8_kind) or real(r4_kind)"
  end select
end function do_time_sum_wrapper

!> Finishes calculations for any reductions that use an average (avg, rms, pow)
!! TODO add mask and any other needed args for adjustment, and pass in the adjusted mask
!! to time_update_done
function diag_reduction_done_wrapper(this, reduction_method, missing_value, has_mask) & !! , has_halo, mask) &
  result(err_msg)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this !< Updated buffer object
  integer, intent(in)                            :: reduction_method !< enumerated reduction type from diag_data
  real(kind=r8_kind), intent(in)                 :: missing_value !< missing_value for masked data points
  logical, intent(in)                            :: has_mask !< indicates if there was a mask used during buffer updates
  character(len=51)                              :: err_msg !< error message to return, blank if sucessful
  !logical, intent(in)                            :: is_subregional !< if subregional output
  !logical, intent(in)                            :: has_halo !< true if halo region is being used
  !logical, optional, intent(in)                  :: mask(:,:,:,:) !< whether a mask variant reduction
  !logical, allocatable                           :: mask_tmp(:,:,:,:)
  !integer :: is, ie, js, je, ks, ke, zs, ze
  !integer :: i, halo_size(4)

  if(.not. allocated(this%buffer)) return

  if(this%weight_sum .eq. 0.0_r8_kind) return

  ! TODO mask adjustment for halos, not needed unless were passing in the mask
  ! if the mask is stil bigger than the buffer, theres a halo region we can leave out
  !if(has_halo .and. present(mask)) then
    !is = lbound(this%buffer,1); ie = ubound(this%buffer,1)
    !js = lbound(this%buffer,2); je = ubound(this%buffer,2)
    !ks = lbound(this%buffer,3); ke = ubound(this%buffer,3)
    !zs = lbound(this%buffer,4); ze = ubound(this%buffer,4)
    !! might be safe to assume these are all the same
    !do i=1, 4
      !halo_size(i) = (SIZE(this%buffer,i) - SIZE(mask,i)) / 2
    !enddo
    !mask_tmp = mask(is+halo_size(1):ie+halo_size(1), js+halo_size(2):je+halo_size(2), ks+halo_size(3):ke+halo_size(3),&
                    !zs+halo_size(4):ze+halo_size(4))
  !endif

  err_msg = ""
  select type(buff => this%buffer)
    type is (real(r8_kind))
      call time_update_done(buff, this%weight_sum, reduction_method, missing_value, has_mask)
    type is (real(r4_kind))
      call time_update_done(buff, this%weight_sum, reduction_method, real(missing_value, r4_kind), has_mask)
  end select
  this%weight_sum = 0.0_r8_kind

end function

!> this leaves out the diurnal index cause its only used for tmp mask allocation
pure function get_buffer_dims(this)
  class(fmsDiagOutputBuffer_type), intent(in) :: this
  integer :: get_buffer_dims(4)
  get_buffer_dims = this%buffer_dims(1:4)
end function

#endif
end module fms_diag_output_buffer_mod
