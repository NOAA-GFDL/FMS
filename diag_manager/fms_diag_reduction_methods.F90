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

!> @defgroup fms_diag_reduction_methods_mod fms_diag_reduction_methods_mod
!> @ingroup diag_manager
!! @brief fms_diag_reduction_methods_mod contains routines that are meant to be used for
!! error checking and setting up to do the reduction methods

!> @file
!> @brief File for @ref fms_diag_reduction_methods_mod

!> @addtogroup fms_diag_reduction_methods_mod
!> @{
module fms_diag_reduction_methods_mod
  use platform_mod, only: r8_kind, r4_kind
  use fms_diag_bbox_mod, only: fmsDiagIbounds_type
  use fms_string_utils_mod, only: string
  use mpp_mod
  implicit none
  private

  public :: check_indices_order, init_mask, set_weight
  public :: do_time_none, do_time_min, do_time_max, do_time_sum_update, time_update_done

  !> @brief Does the time_none reduction method. See include/fms_diag_reduction_methods.inc
  !TODO This needs to be extended to integers
  interface do_time_none
    module procedure do_time_none_r4, do_time_none_r8
  end interface do_time_none

  !> @brief Does the time_min reduction method. See include/fms_diag_reduction_methods.inc
  !TODO This needs to be extended to integers
  interface do_time_min
    module procedure do_time_min_r4, do_time_min_r8
  end interface do_time_min

  !> @brief Does the time_max reduction method. See include/fms_diag_reduction_methods.inc
  !TODO This needs to be extended to integers
  interface do_time_max
    module procedure do_time_max_r4, do_time_max_r8
  end interface do_time_max

  !> @brief Sum update updates the buffer for any reductions that involve summation
  !! (ie. time_sum, avg, rms, pow)
  !!TODO This needs to be extended to integers
  interface do_time_sum_update
    module procedure do_time_sum_update_r4, do_time_sum_update_r8
  end interface

  !> @brief Finishes a reduction that involves an average
  !! (ie. time_avg, rms, pow)
  !! This takes the average at the end of the time step
  interface time_update_done
    module procedure sum_update_done_r4, sum_update_done_r8
  end interface

  contains

  !> @brief Checks improper combinations of is, ie, js, and je.
  !! @return The error message, empty string if no errors were found
  !> @note accept_data works in either one or another of two modes.
  !! 1. Input field is a window (e.g. FMS physics)
  !! 2. Input field includes halo data
  !! It cannot handle a window of data that has halos.
  !! (A field with no windows or halos can be thought of as a special case of either mode.)
  !! The logic for indexing is quite different for these two modes, but is not clearly separated.
  !! If both the beggining and ending indices are present, then field is assumed to have halos.
  !! If only beggining indices are present, then field is assumed to be a window.
  !> @par
  !! There are a number of ways a user could mess up this logic, depending on the combination
  !! of presence/absence of is,ie,js,je. The checks below should catch improper combinations.
  pure function check_indices_order(is_in, ie_in, js_in, je_in) &
  result(error_msg)
    integer, intent(in), optional :: is_in, ie_in, js_in, je_in !< Indices passed to fms_diag_accept_data()
    character(len=128) :: error_msg !< An error message used only for testing purpose!!!

    error_msg = ""
    IF ( PRESENT(ie_in) ) THEN
      IF ( .NOT.PRESENT(is_in) ) THEN
        error_msg = 'ie_in present without is_in'
        return
      END IF
      IF ( PRESENT(js_in) .AND. .NOT.PRESENT(je_in) ) THEN
        error_msg = 'is_in and ie_in present, but js_in present without je_in'
        return
      END IF
    END IF

    IF ( PRESENT(je_in) ) THEN
      IF ( .NOT.PRESENT(js_in) ) THEN
        error_msg = 'je_in present without js_in'
        return
      END IF
      IF ( PRESENT(is_in) .AND. .NOT.PRESENT(ie_in) ) THEN
        error_msg = 'js_in and je_in present, but is_in present without ie_in'
        return
      END IF
    END IF
  end function check_indices_order

  !> @brief Sets the logical mask based on mask or rmask
  !> @return logical mask
  function init_mask(rmask, mask, field) &
  result(oor_mask)
    LOGICAL,  DIMENSION(:,:,:,:), allocatable, INTENT(in) :: mask  !< The location of the mask
    CLASS(*), DIMENSION(:,:,:,:), allocatable, INTENT(in) :: rmask !< The masking values
    CLASS(*), DIMENSION(:,:,:,:),          intent(in) :: field !< Field_data

    logical, allocatable, dimension(:,:,:,:) :: oor_mask !< mask

    ALLOCATE(oor_mask(SIZE(field, 1), SIZE(field, 2), SIZE(field, 3), SIZE(field, 4)))
    oor_mask = .true.

    if (allocated(mask)) then
      oor_mask = mask
    elseif (allocated(rmask)) then
      select type (rmask)
      type is (real(kind=r8_kind))
        WHERE (rmask < 0.5_r8_kind) oor_mask = .FALSE.
      type is (real(kind=r4_kind))
        WHERE (rmask < 0.5_r4_kind) oor_mask = .FALSE.
      end select
    endif

  end function init_mask

  !> @brief Sets the weight based on the weight passed into send_data (1.0_r8_kind if the weight is not passed in)
  !! The weight will be saved as an r8 and converted to r4 as needed
  !! @return weight to use when averaging
  pure function set_weight(weight) &
  result(out_weight)
    CLASS(*), INTENT(in), OPTIONAL :: weight !< The weight use when averaging

    real(kind=r8_kind) :: out_weight

    out_weight = 1.0_r8_kind
    if (present(weight)) then
      select type(weight)
      type is (real(kind=r8_kind))
        out_weight = real(weight, kind = r8_kind)
      type is (real(kind=r4_kind))
        out_Weight = real(weight, kind = r8_kind)
      end select
    endif
  end function set_weight

#include "fms_diag_reduction_methods_r4.fh"
#include "fms_diag_reduction_methods_r8.fh"

end module fms_diag_reduction_methods_mod
!> @}
! close documentation grouping