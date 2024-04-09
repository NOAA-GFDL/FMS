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
use time_manager_mod, only: time_type, operator(==), operator(>=), get_ticks_per_second, get_time, operator(>)
use constants_mod, only: SECONDS_PER_DAY
use mpp_mod, only: mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe
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
  real(r8_kind), allocatable :: weight_sum(:,:,:,:) !< Weight sum as an array
                                                    !! (this will be have a size of 1,1,1,1 when not using variable
                                                    !! masks!)
  integer,  allocatable :: num_elements(:)    !< used in time-averaging
  integer,  allocatable :: axis_ids(:)        !< Axis ids for the buffer
  integer               :: field_id           !< The id of the field the buffer belongs to
  integer               :: yaml_id            !< The id of the yaml id the buffer belongs to
  logical               :: done_with_math     !< .True. if done doing the math
  integer               :: diurnal_sample_size = -1 !< dirunal sample size as read in from the reduction method
                                                    !! ie. diurnal24 = sample size of 24
  integer               :: diurnal_section= -1 !< the diurnal section (ie 5th index) calculated from the current model
                                              !! time and sample size if using a diurnal reduction
  logical, allocatable  :: send_data_called   !< .True. if send_data has been called
  integer               :: unlmited_dimension !< Unlimited dimension index of the last write for this output buffer
  type(time_type)       :: time               !< The last time the data was received
  type(time_type)       :: next_output        !< The next time to output the data

  contains
  procedure :: add_axis_ids
  procedure :: get_axis_ids
  procedure :: set_field_id
  procedure :: get_field_id
  procedure :: set_yaml_id
  procedure :: get_yaml_id
  procedure :: init_buffer_time
  procedure :: set_next_output
  procedure :: update_buffer_time
  procedure :: is_there_data_to_write
  procedure :: is_time_to_finish_reduction
  procedure :: set_send_data_called
  procedure :: is_done_with_math
  procedure :: set_done_with_math
  procedure :: write_buffer
  procedure :: init_buffer_unlim_dim
  procedure :: increase_unlim_dim
  procedure :: get_unlim_dim
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
  procedure :: get_diurnal_sample_size
  procedure :: set_diurnal_sample_size
  procedure :: set_diurnal_section_index
  procedure :: get_remapped_diurnal_data
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
  if (allocated(this%weight_sum))   deallocate(this%weight_sum)
end subroutine flush_buffer

!> Allocates a 5D buffer to given buff_type.
subroutine allocate_buffer(this, buff_type, ndim, buff_sizes, mask_variant, field_name, diurnal_samples)
  class(fmsDiagOutputBuffer_type), intent(inout), target :: this            !< 5D buffer object
  class(*),                        intent(in)            :: buff_type       !< allocates to the type of buff_type
  integer,                         intent(in)            :: ndim            !< Number of dimension
  integer,                         intent(in)            :: buff_sizes(4)   !< dimension buff_sizes
  logical,                         intent(in)            :: mask_variant    !< Mask changes over time
  character(len=*),                intent(in)            :: field_name      !< field name for error output
  integer,                         intent(in)            :: diurnal_samples !< number of diurnal samples

  integer :: n_samples !< number of diurnal samples, defaults to 1

  n_samples = MAX(1, diurnal_samples)
  call this%set_diurnal_sample_size(n_samples)

  this%ndim =ndim
  if(allocated(this%buffer)) call mpp_error(FATAL, "allocate_buffer: buffer already allocated for field:" // &
                                                   field_name)
  select type (buff_type)
    type is (integer(kind=i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                                  & n_samples))
      this%buffer_type = i4
    type is (integer(kind=i8_kind))
      allocate(integer(kind=i8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                                  & n_samples))
      this%buffer_type = i8
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4),  &
                                               & n_samples))
      this%buffer_type = r4
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4), &
                                               & n_samples))
      this%buffer_type = r8
    class default
       call mpp_error("allocate_buffer", &
           "The buff_type value passed to allocate a buffer is not a r8, r4, i8, or i4" // &
           "for field:" // field_name, FATAL)
  end select
  if (mask_variant) then
    allocate(this%weight_sum(buff_sizes(1),buff_sizes(2),buff_sizes(3),buff_sizes(4)))
  else
    allocate(this%weight_sum(1,1,1,1))
  endif
  this%weight_sum = 0.0_r8_kind

  allocate(this%num_elements(n_samples))
  this%num_elements = 0
  this%done_with_math = .false.
  this%send_data_called = .false.
  allocate(this%buffer_dims(5))
  this%buffer_dims(1) = buff_sizes(1)
  this%buffer_dims(2) = buff_sizes(2)
  this%buffer_dims(3) = buff_sizes(3)
  this%buffer_dims(4) = buff_sizes(4)
  this%buffer_dims(5) = n_samples
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
subroutine get_axis_ids(this, res)
  class(fmsDiagOutputBuffer_type), target, intent(inout) :: this        !< Buffer object
  integer, pointer, intent(out) :: res(:)

  if (allocated(this%axis_ids)) then
    res => this%axis_ids
  else
    allocate(res(1))
    res = diag_null
  endif
end subroutine

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
    this%next_output = time
  else
    this%time = get_base_time()
    this%next_output = this%time
  endif
end subroutine init_buffer_time

!> @brief Sets the next output
subroutine set_next_output(this, next_output, next_next_output, is_static)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this             !< Buffer object
  type(time_type),                 intent(in)    :: next_output      !< The current next_output in the file obj
  type(time_type),                 intent(in)    :: next_next_output !< The current next_next_output in the file obj
  logical, optional,               intent(in)    :: is_static        !< .True. if the field is static

  if (present(is_static)) then
    !< If the field is static set the next_output to be equal to time
    !! this should only be used in the init, so next_output will be equal to the the init time
    if (is_static) then
      this%next_output = this%time
      return
    endif
  endif

  !< If the file's next_output is greater than the buffer's next output set
  !! the buffer's next output to the file's next_ouput, otherwise use the file's
  !! next_next_output
  !! This is needed for when file have fields that get data send data sent at different frequencies
  if (next_output > this%next_output) then
    this%next_output = next_output
  else
    this%next_output = next_next_output
  endif
end subroutine set_next_output

!> @brief Update the buffer time if it is a new time
subroutine update_buffer_time(this, time)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  type(time_type),                 intent(in)    :: time        !< Current model time

  if (time > this%time) then
    this%time = time
  endif
end subroutine update_buffer_time

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

!> @brief Get the unlim dimension index of the buffer object
!! @return The unlim dimension index of the buffer object
function get_unlim_dim(this) &
  result(res)
  class(fmsDiagOutputBuffer_type), intent(in) :: this            !< buffer object to write
  integer :: res

  res = this%unlmited_dimension
end function get_unlim_dim

!> @brief Increase the unlim dimension index of the buffer object
subroutine increase_unlim_dim(this)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this            !< buffer object to write

  this%unlmited_dimension = this%unlmited_dimension + 1
end subroutine increase_unlim_dim

!> @brief Init the unlim dimension index of the buffer object to 0
subroutine init_buffer_unlim_dim(this)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this            !< buffer object to write

  this%unlmited_dimension = 0
end subroutine

!> @brief Write the buffer to the file
subroutine write_buffer(this, fms2io_fileobj, unlim_dim_level, is_diurnal)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this            !< buffer object to write
  class(FmsNetcdfFile_t),          intent(in)    :: fms2io_fileobj  !< fileobj to write to
  integer, optional,               intent(in)    :: unlim_dim_level !< unlimited dimension
  logical, optional,               intent(in)    :: is_diurnal !< should be set if using diurnal
                                                       !! reductions so buffer data can be remapped

  select type(fms2io_fileobj)
  type is (FmsNetcdfFile_t)
    call this%write_buffer_wrapper_netcdf(fms2io_fileobj, unlim_dim_level=unlim_dim_level, is_diurnal=is_diurnal)
  type is (FmsNetcdfDomainFile_t)
    call this%write_buffer_wrapper_domain(fms2io_fileobj, unlim_dim_level=unlim_dim_level, is_diurnal=is_diurnal)
  type is (FmsNetcdfUnstructuredDomainFile_t)
    call this%write_buffer_wrapper_u(fms2io_fileobj, unlim_dim_level=unlim_dim_level, is_diurnal=is_diurnal)
  class default
    call mpp_error(FATAL, "The file "//trim(fms2io_fileobj%path)//" is not one of the accepted types"//&
      " only FmsNetcdfFile_t, FmsNetcdfDomainFile_t, and FmsNetcdfUnstructuredDomainFile_t are accepted.")
  end select

  call this%initialize_buffer(diag_yaml%diag_fields(this%yaml_id)%get_var_reduction(), &
    diag_yaml%diag_fields(this%yaml_id)%get_var_outname())
  !TODO Set the counters back to 0
end subroutine write_buffer

!> @brief Write the buffer to the FmsNetcdfFile_t fms2io_fileobj
subroutine write_buffer_wrapper_netcdf(this, fms2io_fileobj, unlim_dim_level, is_diurnal)
  class(fmsDiagOutputBuffer_type),  intent(in) :: this            !< buffer object to write
  type(FmsNetcdfFile_t),            intent(in) :: fms2io_fileobj  !< fileobj to write to
  integer, optional,                intent(in) :: unlim_dim_level !< unlimited dimension
  logical, optional,                intent(in) :: is_diurnal !< should be set if using diurnal
                                                       !! reductions so buffer data can be remapped
  character(len=:), allocatable :: varname !< name of the variable
  logical :: using_diurnal !< local copy of is_diurnal if present
  class(*), allocatable                        :: buff_ptr(:,:,:,:,:) !< pointer for buffer to write

  using_diurnal = .false.
  if( present(is_diurnal) ) using_diurnal = is_diurnal
  if( using_diurnal ) then
    call this%get_remapped_diurnal_data(buff_ptr)
  else
    buff_ptr = this%buffer
  endif

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select case(this%ndim)
  case (0)
    call write_data(fms2io_fileobj, varname, buff_ptr(1,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (1)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (2)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,1,1,1), unlim_dim_level=unlim_dim_level)
  case (3)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,1,1), unlim_dim_level=unlim_dim_level)
  case (4)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,:,1), unlim_dim_level=unlim_dim_level)
  case (5)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,:,:), unlim_dim_level=unlim_dim_level)
  end select
end subroutine write_buffer_wrapper_netcdf

!> @brief Write the buffer to the FmsNetcdfDomainFile_t fms2io_fileobj
subroutine write_buffer_wrapper_domain(this, fms2io_fileobj, unlim_dim_level, is_diurnal)
  class(fmsDiagOutputBuffer_type),    intent(in) :: this            !< buffer object to write
  type(FmsNetcdfDomainFile_t),        intent(in) :: fms2io_fileobj  !< fileobj to write to
  integer, optional,                  intent(in) :: unlim_dim_level !< unlimited dimension
  logical, optional,                  intent(in) :: is_diurnal !< should be set if using diurnal
                                                       !! reductions so buffer data can be remapped

  character(len=:), allocatable :: varname !< name of the variable
  logical :: using_diurnal !< local copy of is_diurnal if present
  class(*), allocatable                        :: buff_ptr(:,:,:,:,:) !< pointer to buffer to write

  using_diurnal = .false.
  if( present(is_diurnal) ) using_diurnal = is_diurnal
  if( using_diurnal ) then
    call this%get_remapped_diurnal_data(buff_ptr)
  else
    buff_ptr = this%buffer
  endif

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select case(this%ndim)
  case (0)
    call write_data(fms2io_fileobj, varname, buff_ptr(1,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (1)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (2)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,1,1,1), unlim_dim_level=unlim_dim_level)
  case (3)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,1,1), unlim_dim_level=unlim_dim_level)
  case (4)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,:,1), unlim_dim_level=unlim_dim_level)
  case (5)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,:,:), unlim_dim_level=unlim_dim_level)
  end select
end subroutine write_buffer_wrapper_domain

!> @brief Write the buffer to the FmsNetcdfUnstructuredDomainFile_t fms2io_fileobj
subroutine write_buffer_wrapper_u(this, fms2io_fileobj, unlim_dim_level, is_diurnal)
  class(fmsDiagOutputBuffer_type),                 intent(in) :: this            !< buffer object to write
  type(FmsNetcdfUnstructuredDomainFile_t),         intent(in) :: fms2io_fileobj  !< fileobj to write to
  integer, optional,                               intent(in) :: unlim_dim_level !< unlimited dimension
  logical, optional,                  intent(in) :: is_diurnal !< should be set if using diurnal
                                                       !! reductions so buffer data can be remapped

  character(len=:), allocatable :: varname !< name of the variable
  logical :: using_diurnal !< local copy of is_diurnal if present
  class(*), allocatable                        :: buff_ptr(:,:,:,:,:) !< pointer for buffer to write

  using_diurnal = .false.
  if( present(is_diurnal) ) using_diurnal = is_diurnal
  if( using_diurnal ) then
    call this%get_remapped_diurnal_data(buff_ptr)
  else
    buff_ptr = this%buffer
  endif

  varname = diag_yaml%diag_fields(this%yaml_id)%get_var_outname()
  select case(this%ndim)
  case (0)
    call write_data(fms2io_fileobj, varname, buff_ptr(1,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (1)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,1,1,1,1), unlim_dim_level=unlim_dim_level)
  case (2)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,1,1,1), unlim_dim_level=unlim_dim_level)
  case (3)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,1,1), unlim_dim_level=unlim_dim_level)
  case (4)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,:,1), unlim_dim_level=unlim_dim_level)
  case (5)
    call write_data(fms2io_fileobj, varname, buff_ptr(:,:,:,:,:), unlim_dim_level=unlim_dim_level)
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
function do_time_sum_wrapper(this, field_data, mask, is_masked, mask_variant, bounds_in, bounds_out, missing_value, &
                             has_missing_value, pow_value) &
  result(err_msg)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this                !< buffer object to write
  class(*),                        intent(in)    :: field_data(:,:,:,:) !< Buffer data for current time
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_in           !< Indicies for the buffer passed in
  type(fmsDiagIbounds_type),       intent(in)    :: bounds_out          !< Indicies for the output buffer
  logical,                         intent(in)    :: mask(:,:,:,:)       !< Mask for the field
  logical,                         intent(in)    :: is_masked           !< .True. if the field has a mask
  logical,                         intent(in)    :: mask_variant        !< .True. if the mask changes over time
  real(kind=r8_kind),              intent(in)    :: missing_value       !< Missing_value for data points that are masked
  logical,                         intent(in)    :: has_missing_value   !< .True. if the field was registered with
                                                                        !! a missing value
  integer, optional,               intent(in)    :: pow_value           !< power value, will calculate field_data^pow
                                                                        !! before adding to buffer should only be
                                                                        !! present if using pow reduction method
  character(len=150) :: err_msg

  !TODO This will be expanded for integers
  err_msg = ""
  select type (output_buffer => this%buffer)
    type is (real(kind=r8_kind))
      select type (field_data)
      type is (real(kind=r8_kind))
        if (.not. is_masked) then
          if (any(field_data .eq. missing_value)) &
            err_msg = "You cannot pass data with missing values without masking them!"
        endif
        call do_time_sum_update(output_buffer, this%weight_sum, field_data, mask, is_masked, mask_variant, &
                                bounds_in, bounds_out, missing_value, this%diurnal_section, &
                                pow=pow_value)
      class default
        err_msg="do_time_sum_wrapper::the output buffer and the buffer send in are not of the same type (r8_kind)"
      end select
    type is (real(kind=r4_kind))
      select type (field_data)
      type is (real(kind=r4_kind))
        if (.not. is_masked) then
          if (any(field_data .eq. missing_value)) &
            err_msg = "You cannot pass data with missing values without masking them!"
        endif
        call do_time_sum_update(output_buffer, this%weight_sum, field_data, mask, is_masked, mask_variant, &
                                bounds_in, bounds_out, real(missing_value, kind=r4_kind), &
                                this%diurnal_section, pow=pow_value)
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
function diag_reduction_done_wrapper(this, reduction_method, missing_value, has_mask, mask_variant) &
  result(err_msg)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this !< Updated buffer object
  integer, intent(in)                            :: reduction_method !< enumerated reduction type from diag_data
  real(kind=r8_kind), intent(in)                 :: missing_value !< missing_value for masked data points
  logical, intent(in)                            :: has_mask !< indicates if there was a mask used during buffer updates
  logical, intent(in)                            :: mask_variant !< Indicates if the mask changes over time
  character(len=51)                              :: err_msg !< error message to return, blank if sucessful

  if(.not. allocated(this%buffer)) return

  err_msg = ""
  select type(buff => this%buffer)
    type is (real(r8_kind))
      call time_update_done(buff, this%weight_sum, reduction_method, missing_value, has_mask, mask_variant, &
        this%diurnal_sample_size)
    type is (real(r4_kind))
      call time_update_done(buff, this%weight_sum, reduction_method, real(missing_value, r4_kind), has_mask, &
                            mask_variant, this%diurnal_sample_size)
  end select
  this%weight_sum = 0.0_r8_kind

end function

!> this leaves out the diurnal index cause its only used for tmp mask allocation
pure function get_buffer_dims(this)
  class(fmsDiagOutputBuffer_type), intent(in) :: this !< buffer object to get from
  integer :: get_buffer_dims(4)
  get_buffer_dims = this%buffer_dims(1:4)
end function

!> Get diurnal sample size (amount of diurnal sections)
pure integer function get_diurnal_sample_size(this)
  class(fmsDiagOutputBuffer_type), intent(in) :: this !< buffer object to get from
  get_diurnal_sample_size = this%diurnal_sample_size
end function get_diurnal_sample_size

!> Set diurnal sample size (amount of diurnal sections)
subroutine set_diurnal_sample_size(this, sample_size)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this !< buffer object to set sample size for
  integer, intent(in)                            :: sample_size !< sample size to used to split daily
                                                               !! data into given amount of sections
  this%diurnal_sample_size = sample_size
end subroutine set_diurnal_sample_size

!> Set diurnal section index based off the current time and previously set diurnal_samplesize
!! Calculates which diurnal section of daily data the current time is in
subroutine set_diurnal_section_index(this, time)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this !< buffer object to set diurnal index for
  type(time_type), intent(in)                     :: time !< current model time
  integer :: seconds, days, ticks

  if(this%diurnal_sample_size .lt. 0) call mpp_error(FATAL, "set_diurnal_section_index::"// &
    " diurnal sample size must be set before trying to set diurnal index for send_data")

  call get_time(time,seconds,days,ticks) ! get current date
  ! calculates which diurnal section current time is in for a given amount of diurnal sections(<24)
  this%diurnal_section = floor( (seconds+real(ticks)/get_ticks_per_second()) &
                       & * this%diurnal_sample_size/SECONDS_PER_DAY) + 1
end subroutine set_diurnal_section_index

!> Remaps the output buffer array when using the diurnal reduction
!! moves the diurnal index to the left-most unused dimension for the io
subroutine get_remapped_diurnal_data(this, res)
  class(fmsDiagOutputBuffer_type), intent(in) :: this !< output buffer object
  class(*), intent(out), allocatable :: res(:,:,:,:,:) !< resulting remapped data
  integer :: last_dim !< last dimension thats used
  integer :: ie, je, ke, ze, de !< ending indices for the new array
  integer(i4_kind) :: buff_size(5)!< sizes for allocated buffer

  ! last dim is number of dimensions - 1 for diurnal axis
  last_dim = this%ndim - 1
  ! get the bounds of the remapped output array based on # of dims
  ke = 1; ze = 1; de = 1
  select case(last_dim)
    case (1)
      ie = this%buffer_dims(1); je = this%buffer_dims(5)
    case (2)
      ie = this%buffer_dims(1); je = this%buffer_dims(2)
      ke = this%buffer_dims(5)
    case (3)
      ie = this%buffer_dims(1); je = this%buffer_dims(2)
      ke = this%buffer_dims(3); ze = this%buffer_dims(5)
    case (4)
      ! no need to remap if 4d
      res = this%buffer
      return
  end select

  select type(buff => this%buffer)
    type is (real(r8_kind))
      allocate(real(r8_kind) :: res(1:ie, 1:je, 1:ke, 1:ze, 1:de))
      select type(res)
        type is (real(r8_kind))
          res(1:ie, 1:je, 1:ke, 1:ze, 1:de) = reshape(buff, SHAPE(res))
      end select
    type is (real(r4_kind))
      allocate(real(r4_kind) :: res(1:ie, 1:je, 1:ke, 1:ze, 1:de))
      select type(res)
        type is (real(r4_kind))
          res(1:ie, 1:je, 1:ke, 1:ze, 1:de) = reshape(buff, SHAPE(res))
      end select
    type is (integer(i8_kind))
      allocate(integer(i8_kind) :: res(1:ie, 1:je, 1:ke, 1:ze, 1:de))
      select type(res)
        type is (integer(i8_kind))
          res(1:ie, 1:je, 1:ke, 1:ze, 1:de) = reshape(buff, SHAPE(res))
      end select
    type is (integer(i4_kind))
      allocate(integer(i4_kind) :: res(1:ie, 1:je, 1:ke, 1:ze, 1:de))
      select type(res)
        type is (integer(i4_kind))
          res(1:ie, 1:je, 1:ke, 1:ze, 1:de) = reshape(buff, SHAPE(res))
      end select
  end select

end subroutine get_remapped_diurnal_data

!> @brief Determine if there is any data to write (i.e send_data has been called)
!! @return .true. if there is data to write
function is_there_data_to_write(this) &
  result(res)
  class(fmsDiagOutputBuffer_type), intent(in) :: this        !< Buffer object

  logical :: res

  if (allocated(this%send_data_called)) then
    res = this%send_data_called
  else
    res = .false.
  endif
end function

!> @brief Determine if it is time to finish the reduction method
!! @return .true. if it is time to finish the reduction method
function is_time_to_finish_reduction(this, end_time) &
  result(res)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object
  type(time_type), optional,       intent(in)    :: end_time    !< The time at the end of the run

  logical :: res

  res = .false.
  if (this%time >= this%next_output) res = .true.

  if (present(end_time)) then
    if (end_time >= this%next_output) res = .true.
  endif
end function is_time_to_finish_reduction

!> @brief Sets send_data_called to .true.
subroutine set_send_data_called(this)
  class(fmsDiagOutputBuffer_type), intent(inout) :: this        !< Buffer object

  this%send_data_called = .true.
end subroutine set_send_data_called
#endif
end module fms_diag_output_buffer_mod
