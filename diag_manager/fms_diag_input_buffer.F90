!***********************************************************************
!*                   GNU Lesser General Public License
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
!> @defgroup fms_diag_input_buffer_mod fms_diag_input_buffer_mod
!> @ingroup diag_manager
!! @brief
!> @addtogroup fms_diag_input_buffer_mod
!> @{
module fms_diag_input_buffer_mod
#ifdef use_yaml
  use platform_mod,             only: r8_kind, r4_kind, i4_kind, i8_kind
  use fms_diag_axis_object_mod, only: fmsDiagAxisContainer_type, fmsDiagFullAxis_type
  use time_manager_mod,         only: time_type
  use mpp_mod,                  only: mpp_error, FATAL
  implicit NONE
  private

  !> @brief Appends the input_data_buffer and the mask (only when the mask is set to .True.)
  interface append_data_buffer
    module procedure append_data_buffer_r4, append_data_buffer_r8
  end interface

  !> @brief Sums the data in the input_data_buffer
  interface sum_data_buffer
    module procedure sum_data_buffer_r4, sum_data_buffer_r8
  end interface

  !> @brief Type to hold the information needed for the input buffer
  !! This is used when set_math_needs_to_be_done = .true. (i.e calling send_data
  !! from an openmp region with multiple threads)
  type fmsDiagInputBuffer_t
    logical                        :: initialized     !< .True. if the input buffer has been initialized
    class(*),          allocatable :: buffer(:,:,:,:) !< Input data passed in send_data
    integer,           allocatable :: counter(:,:,:,:)!< Number of send_data calls for each point
    real(kind=r8_kind)             :: weight          !< Weight passed in send_data
    type(time_type)                :: send_data_time  !< The time send data was called last

    contains
    procedure :: get_buffer
    procedure :: get_weight
    procedure :: allocate_input_buffer_object
    procedure :: init_input_buffer_object
    procedure :: set_input_buffer_object
    procedure :: update_input_buffer_object
    procedure :: prepare_input_buffer_object
    procedure :: set_send_data_time
    procedure :: get_send_data_time
    procedure :: is_initialized
  end type fmsDiagInputBuffer_t

  public :: fmsDiagInputBuffer_t

  contains

  !> @brief Get the buffer from the input buffer object
  !! @return a pointer to the buffer
  function get_buffer(this) &
    result(buffer)
    class(fmsDiagInputBuffer_t), target, intent(in) :: this !< input buffer object
    class(*), pointer :: buffer(:,:,:,:)

    buffer => this%buffer
  end function get_buffer


  !> @brief Get the weight from the input buffer object
  !! @return a pointer to the weight
  function get_weight(this) &
    result(weight)
    class(fmsDiagInputBuffer_t), target, intent(in) :: this !< input buffer object
    real(kind=r8_kind), pointer :: weight

    weight => this%weight
  end function get_weight

  !> @brief Initiliazes an input data buffer
  !! @return Error message if something went wrong
  function allocate_input_buffer_object(this, input_data, axis_ids, diag_axis) &
    result(err_msg)
    class(fmsDiagInputBuffer_t),         intent(out)   :: this                !< input buffer object
    class(*),                            intent(in)    :: input_data(:,:,:,:) !< input data
    integer, target,                     intent(in)    :: axis_ids(:)         !< axis ids for the field
    class(fmsDiagAxisContainer_type),    intent(in)    :: diag_axis(:)        !< Array of diag_axis
    character(len=128) :: err_msg

    integer            :: naxes         !< The number of axes in the field
    integer, parameter :: ndims = 4     !< Number of dimensions
    integer            :: length(ndims) !< The length of an axis
    integer            :: a             !< For looping through axes
    integer, pointer   :: axis_id       !< The axis ID

    err_msg = ""

    !! Use the axis to get the size
    !> Initialize the axis lengths to 1.  Any dimension that does not have an axis will have a length
    !! of 1.
    length = 1
    naxes = size(axis_ids)
    axis_loop: do a = 1,naxes
      axis_id => axis_ids(a)
      select type (axis => diag_axis(axis_id)%axis)
      type is (fmsDiagFullAxis_type)
        length(a) = axis%axis_length()
      end select
    enddo axis_loop

    select type (input_data)
    type is (real(r4_kind))
      allocate(real(kind=r4_kind) :: this%buffer(length(1), length(2), length(3), length(4)))
      this%buffer = 0.0_r4_kind
    type is (real(r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(length(1), length(2), length(3), length(4)))
      this%buffer = 0.0_r8_kind
    type is (integer(i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(length(1), length(2), length(3), length(4)))
      this%buffer = 0_i4_kind
    type is (integer(i8_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(length(1), length(2), length(3), length(4)))
      this%buffer = 0_i8_kind
    class default
      err_msg = "The data input is not one of the supported types. &
                &Only r4, r8, i4, and i8 types are supported."
    end select

    this%weight = 1.0_r8_kind
    this%initialized = .true.
    allocate(this%counter(length(1), length(2), length(3), length(4)))
    this%counter = 0
  end function allocate_input_buffer_object

  !> @brief Initiliazes an input data buffer and the counter
  subroutine init_input_buffer_object(this)
    class(fmsDiagInputBuffer_t), intent(inout) :: this                !< input buffer object

    select type(buffer=>this%buffer)
    type is (real(kind=r8_kind))
      buffer = 0.0_r8_kind
    type is (real(kind=r4_kind))
      buffer = 0.0_r4_kind
    end select
    this%counter = 0
  end subroutine init_input_buffer_object

  !> @brief Sets the time send data was called last
  subroutine set_send_data_time(this, time)
    class(fmsDiagInputBuffer_t), intent(inout) :: this                !< input buffer object
    type(time_type),             intent(in)    :: time                !< The time send data was called

    this%send_data_time = time
  end subroutine set_send_data_time

  !> @brief Get the time send data was called last
  !! @result the time send data was called last
  function get_send_data_time(this) &
    result(rslt)
    class(fmsDiagInputBuffer_t), intent(in) :: this                !< input buffer object
    type(time_type) :: rslt

    rslt = this%send_data_time
  end function get_send_data_time

  !> @brief Updates the input data buffer object for the current send_data call
  !! @return Error message (if an error occurs)
  function update_input_buffer_object(this, input_data, is, js, ks, ie, je, ke, mask_in, mask_out, &
                                      mask_variant, var_is_masked) &
    result(err_msg)

    class(fmsDiagInputBuffer_t), intent(inout) :: this                !< input buffer object
    class(*),                    intent(in)    :: input_data(:,:,:,:) !< Field data
    integer,                     intent(in)    :: is, js, ks          !< Starting index for each of the dimension
    integer,                     intent(in)    :: ie, je, ke          !< Ending index for each of the dimensions
    logical,                     intent(in)    :: mask_in(:,:,:,:)
    logical,                     intent(inout) :: mask_out(:,:,:,:)
    logical,                     intent(in)    :: mask_variant
    logical,                     intent(in)    :: var_is_masked

    character(len=128) :: err_msg

    if (mask_variant) then
      err_msg = append_data_buffer_wrapper(mask_out(is:ie,js:je,ks:ke,:), mask_in, &
                        this%buffer(is:ie,js:je,ks:ke,:), input_data)
    else
      mask_out(is:ie,js:je,ks:ke,:) = mask_in
      err_msg = sum_data_buffer_wrapper(mask_in, this%buffer(is:ie,js:je,ks:ke,:), input_data, &
                                   this%counter(is:ie,js:je,ks:ke,:), &
                                   var_is_masked)
    endif

  end function update_input_buffer_object

  !> @brief Prepare the input data buffer to do the reduction methods (i.e divide by the number of times
  !! send data has been called)
  subroutine prepare_input_buffer_object(this, field_info)
    class(fmsDiagInputBuffer_t), intent(inout) :: this                !< input buffer object
    character(len=*),            intent(in)    :: field_info          !< Field info to append to error message

    select type (input_data => this%buffer)
    type is (real(kind=r4_kind))
      input_data = input_data / this%counter(1,1,1,1)
    type is (real(kind=r8_kind))
      input_data = input_data / this%counter(1,1,1,1)
    class default
      call mpp_error(FATAL, "prepare_input_buffer_object::"//trim(field_info)//&
                            " has only been implemented for real variables. Contact developers.")
    end select
  end subroutine prepare_input_buffer_object

  !> @brief Sums the data in the input_data_buffer
  !! @return Error message (if an error occurs)
  function sum_data_buffer_wrapper(mask, data_out, data_in, counter, var_is_masked) &
    result(err_msg)

    logical,  intent(in)    :: mask(:,:,:,:)     !< Mask passed into send_data
    class(*), intent(inout) :: data_out(:,:,:,:) !< Data currently saved in the input_data_buffer
    class(*), intent(in)    :: data_in(:,:,:,:)  !< Data passed into send_data
    integer,  intent(inout) :: counter(:,:,:,:)  !< Number of times data has been summed
    logical,  intent(in)    :: var_is_masked     !< .True. if the variable is masked

    character(len=128) :: err_msg

    err_msg = ""
    select type(data_out)
    type is (real(kind=r8_kind))
      select type (data_in)
      type is (real(kind=r8_kind))
        call sum_data_buffer(mask, data_out, data_in, counter, var_is_masked)
      end select
    type is (real(kind=r4_kind))
      select type (data_in)
      type is (real(kind=r4_kind))
        call sum_data_buffer(mask, data_out, data_in, counter, var_is_masked)
      end select
    class default
      err_msg = "sum_data_buffer_wrapper:: has only been implemented for real. Contact developers"
    end select
  end function sum_data_buffer_wrapper

  !> @brief Appends the input_data_buffer and the mask (only when the mask is set to .True.)
  !! @return Error message (if an error occurs)
  function append_data_buffer_wrapper(mask_out, mask_in, data_out, data_in) &
    result(err_msg)
    logical,  intent(inout) :: mask_out(:,:,:,:) !< Mask currently in the input_data_buffer
    logical,  intent(in)    :: mask_in(:,:,:,:)  !< Mask passed in to send_data
    class(*), intent(inout) :: data_out(:,:,:,:) !< Data currently in the input_data_buffer
    class(*), intent(in)    :: data_in(:,:,:,:)  !< Data passed in to send_data

    character(len=128) :: err_msg

    err_msg = ""
    select type(data_out)
    type is (real(kind=r8_kind))
      select type (data_in)
      type is (real(kind=r8_kind))
        call append_data_buffer(mask_out, mask_in, data_out, data_in)
      end select
    type is (real(kind=r4_kind))
      select type (data_in)
      type is (real(kind=r4_kind))
        call append_data_buffer(mask_out, mask_in, data_out, data_in)
      end select
    class default
      err_msg = "append_data_buffer:: has only been implemented for real. Contact developers"
    end select
  end function append_data_buffer_wrapper

  !> @brief Sets the members of the input buffer object
  !! @return Error message if something went wrong
  function set_input_buffer_object(this, input_data, weight, is, js, ks, ie, je, ke) &
    result(err_msg)

    class(fmsDiagInputBuffer_t), intent(inout) :: this                !< input buffer object
    class(*),                    intent(in)    :: input_data(:,:,:,:) !< Field data
    real(kind=r8_kind),          intent(in)    :: weight              !< Weight for the field
    integer,                     intent(in)    :: is, js, ks          !< Starting index for each of the dimension
    integer,                     intent(in)    :: ie, je, ke          !< Ending index for each of the dimensions

    character(len=128) :: err_msg
    err_msg = ""

    if (.not. this%initialized) then
      err_msg = "The data buffer was never initiliazed. This shouldn't happen."
      return
    endif

    this%weight = weight

    select type (input_data)
    type is (real(kind=r4_kind))
      select type (db => this%buffer)
      type is (real(kind=r4_kind))
        db(is:ie, js:je, ks:ke, :) = input_data
      class default
        err_msg = "The data buffer was not allocated to the correct type (r4_kind). This shouldn't happen"
        return
      end select
    type is (real(kind=r8_kind))
      select type (db => this%buffer)
      type is (real(kind=r8_kind))
        db(is:ie, js:je, ks:ke, :) = input_data
      class default
        err_msg = "The data buffer was not allocated to the correct type (r8_kind). This shouldn't happen"
        return
      end select
    type is (integer(kind=i4_kind))
      select type (db => this%buffer)
      type is (integer(kind=i4_kind))
        db(is:ie, js:je, ks:ke, :) = input_data
      class default
        err_msg = "The data buffer was not allocated to the correct type (i4_kind). This shouldn't happen"
        return
      end select
    type is (integer(kind=i8_kind))
      select type (db => this%buffer)
      type is (integer(kind=i8_kind))
        db(is:ie, js:je, ks:ke, :) = input_data
      class default
        err_msg = "The data buffer was not allocated to the correct type (i8_kind). This shouldn't happen"
        return
      end select
    end select
  end function set_input_buffer_object

  !> @brief Determine if an input buffer is initialized
  !! @return .true. if the input buffer is initialized
  pure logical function is_initialized(this)
    class(fmsDiagInputBuffer_t), intent(in) :: this !< input buffer object

    is_initialized = .false.
    if (this%initialized) then
      is_initialized = .true.
    else
      if (allocated(this%buffer)) is_initialized = .true.
    endif
  end function is_initialized

#include "fms_diag_input_buffer_r4.fh"
#include "fms_diag_input_buffer_r8.fh"

#endif
end module fms_diag_input_buffer_mod
!> @}
