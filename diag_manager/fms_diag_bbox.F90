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

!> @defgroup fms_diag_bbox_mod fms_diag_bbox_mod
!> @ingroup diag_manager
!> @brief fms_diag_bbox_mod defines classes encapsulating bounding boxes
!!   and interval bounds.
!!
!> @author Miguel Zuniga
!!
!> @file
!> @brief File for @ref fms_diag_bbox_mod
!> @addtogroup fms_diag_bbox_mod
!> @{
MODULE fms_diag_bbox_mod

   USE fms_mod, ONLY: error_mesg, FATAL

   implicit none

!> @brief Data structure holding a 3D bounding box. It is commonlyused to
!! represent the interval bounds or limits of a 3D sub-array such as the
!! array index bounds of the spatial component a diag_manager field output
!! buffer array.
   TYPE, public :: fmsDiagIbounds_type
      PRIVATE
      INTEGER :: imin !< Lower i bound.
      INTEGER :: imax !< Upper i bound.
      INTEGER :: jmin !< Lower j bound.
      INTEGER :: jmax !< Upper j bound.
      INTEGER :: kmin !< Lower k bound.
      INTEGER :: kmax !< Upper k bound.
   contains
      procedure :: reset => reset_bounds
      procedure :: reset_bounds_from_array_4D
      procedure :: reset_bounds_from_array_5D
      procedure :: update_bounds
      procedure :: get_imin
      procedure :: get_imax
      procedure :: get_jmin
      procedure :: get_jmax
      procedure :: get_kmin
      procedure :: get_kmax
   END TYPE fmsDiagIbounds_type

CONTAINS

   !> @brief Gets imin of fmsDiagIbounds_type
   !! @return copy of integer member imin
   pure integer function get_imin (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%imin
   end function get_imin

   !> @brief Gets imax of fmsDiagIbounds_type
   !! @return copy of integer member imax
   pure integer function get_imax (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%imax
   end function get_imax

   !> @brief Gets jmin of fmsDiagIbounds_type
   !! @return copy of integer member jmin
   pure integer function get_jmin (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%jmin
   end function get_jmin

   !> @brief Gets jmax of fmsDiagIbounds_type
   !! @return copy of integer member jmax
   pure integer function get_jmax (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%jmax
   end function get_jmax


   !> @brief Gets kmin of fmsDiagIbounds_type
   !! @return copy of integer member kmin
   pure integer function get_kmin (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%kmin
   end function get_kmin

   !> @brief Gets kmax of fmsDiagIbounds_type
   !! @return copy of integer member kmax
   pure integer function get_kmax (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%kmax
   end function get_kmax

   !> @brief Reset the instance bounding lower and upper bounds to lower_val and upper_val, respectively.
   SUBROUTINE reset_bounds (this, lower_val, upper_val)
      class (fmsDiagIbounds_type), target, intent(inout) :: this   !< ibounds instance
      integer, intent(in) :: lower_val  !< value for the lower bounds in each dimension
      integer, intent(in) :: upper_val  !< value for the upper bounds in each dimension
      this%imin = lower_val
      this%jmin = lower_val
      this%kmin = lower_val
      this%imax = upper_val
      this%jmax = upper_val
      this%kmax = upper_val
   END SUBROUTINE reset_bounds

   !> @brief Update the the first three (normally  x, y, and z)  min and
   !! max boundaries (array indices) of the instance bounding box
   !! the six specified bounds values.
   SUBROUTINE update_bounds(this, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k)
      CLASS  (fmsDiagIbounds_type), intent(inout) :: this !<The bounding box of the output field buffer inindex space.
      INTEGER, INTENT(in) :: lower_i !< Lower i bound.
      INTEGER, INTENT(in) :: upper_i !< Upper i bound.
      INTEGER, INTENT(in) :: lower_j !< Lower j bound.
      INTEGER, INTENT(in) :: upper_j !< Upper j bound.
      INTEGER, INTENT(in) :: lower_k !< Lower k bound.
      INTEGER, INTENT(in) :: upper_k !< Upper k bound.
      this%imin = MIN(this%imin, lower_i)
      this%imax = MAX(this%imax, upper_i)
      this%jmin = MIN(this%jmin, lower_j)
      this%jmax = MAX(this%jmax, upper_j)
      this%kmin = MIN(this%kmin, lower_k)
      this%kmax = MAX(this%kmax, upper_k)
   END SUBROUTINE update_bounds

   !> @brief Reset the instance bounding box with the bounds determined from the
   !! first three dimensions of the 5D "array" argument
   SUBROUTINE reset_bounds_from_array_4D(this, array)
      CLASS (fmsDiagIbounds_type), INTENT(inout) :: this !< The instance of the bounding box.
      REAL, INTENT( in), DIMENSION(:,:,:,:) :: array !< The 4D input array.
      this%imin = LBOUND(array,1)
      this%imax = UBOUND(array,1)
      this%jmin = LBOUND(array,2)
      this%jmax = UBOUND(array,2)
      this%kmin = LBOUND(array,3)
      this%kmax = UBOUND(array,3)
   END SUBROUTINE  reset_bounds_from_array_4D

   !> @brief Reset the instance bounding box with the bounds determined from the
   !! first three dimensions of the 5D "array" argument
   SUBROUTINE reset_bounds_from_array_5D(this, array)
      CLASS (fmsDiagIbounds_type), INTENT(inout) :: this !< The instance of the bounding box.
      CLASS(*), INTENT( in), DIMENSION(:,:,:,:,:) :: array !< The 5D input array.
      this%imin = LBOUND(array,1)
      this%imax = UBOUND(array,1)
      this%jmin = LBOUND(array,2)
      this%jmax = UBOUND(array,2)
      this%kmin = LBOUND(array,3)
      this%kmax = UBOUND(array,3)
   END SUBROUTINE  reset_bounds_from_array_5D

  END MODULE fms_diag_bbox_mod
  !> @}
  ! close documentation grouping
