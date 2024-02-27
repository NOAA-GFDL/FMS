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
  implicit NONE
  private

  !> @brief Type to hold the information needed for the input buffer
  !! This is used when set_math_needs_to_be_done = .true. (i.e calling send_data
  !! from an openmp region with multiple threads)
  type fmsDiagInputBuffer_t
    logical                        :: initialized     !< .True. if the input buffer has been initialized
    class(*),          allocatable :: buffer(:,:,:,:) !< Input data passed in send_data
    real(kind=r8_kind)             :: weight          !< Weight passed in send_data
    type(time_type)                :: send_data_time  !< The time send data was called

    contains
    procedure :: get_buffer
    procedure :: get_weight
    procedure :: init => init_input_buffer_object
    procedure :: set_input_buffer_object
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
  function init_input_buffer_object(this, input_data, axis_ids, diag_axis) &
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
    type is (real(r8_kind))
      allocate(real(kind=r8_kind) :: this%buffer(length(1), length(2), length(3), length(4)))
    type is (integer(i4_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(length(1), length(2), length(3), length(4)))
    type is (integer(i8_kind))
      allocate(integer(kind=i4_kind) :: this%buffer(length(1), length(2), length(3), length(4)))
    class default
      err_msg = "The data input is not one of the supported types."&
        "Only r4, r8, i4, and i8 types are supported."
    end select

    this%weight = 1.0_r8_kind
    this%initialized = .true.
  end function init_input_buffer_object

  !> @brief Sets the send data time
  subroutine set_send_data_time(this, time)
    class(fmsDiagInputBuffer_t), intent(inout) :: this                !< input buffer object
    type(time_type),             intent(in)    :: time                !< The time send data was called

    this%send_data_time = time
  end subroutine set_send_data_time

  !> @brief Get the send data time
  !! @result Send data time
  function get_send_data_time(this) &
    result(rslt)
    class(fmsDiagInputBuffer_t), intent(in) :: this                !< input buffer object
    type(time_type) :: rslt

    rslt = this%send_data_time
  end function get_send_data_time

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
#endif
end module fms_diag_input_buffer_mod
!> @}
