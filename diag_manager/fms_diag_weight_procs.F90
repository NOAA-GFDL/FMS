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

!> @defgroup fms_diag_weight_procs_mod fms_diag_weight_procs_mod
!> @ingroup diag_manager
!> @brief fms_diag_weight_procs_mod defines procedures for
!! weighting the output fields.
!!
!> @author Miguel Zuniga
!!
!! <TT>fms_diag_weight_procs_mod</TT>  defines procedures, interfaces and
!! used in weighting output fields. They are useful in setting a pointer
!! to the weighting function and (thereby) moving logic of selecting
!! the function from an inner loop to a (more) outer loop. See, for example,
!! their use in the send_data_3d function.
!!
!> @file
!> @brief File for @ref fms_diag_weight_procs_mod
!> @addtogroup fms_diag_weight_procs_mod
!> @{
MODULE fms_diag_weight_procs_mod
   IMPLICIT NONE

   !> Interface for the scalar field wheighting functions.
   ABSTRACT INTERFACE
      FUNCTION weight_the_field ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL :: weight_the_field
      END FUNCTION weight_the_field
   END INTERFACE

   !> Interface for the 3D field wheighting functions.
   ABSTRACT INTERFACE
      FUNCTION weight_the_field_3d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:,:,:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d
      END FUNCTION weight_the_field_3d
   END INTERFACE

   !> Interface for the 1D field wheighting functions.
   ABSTRACT INTERFACE
      FUNCTION weight_the_field_1d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d
      END FUNCTION weight_the_field_1d
   END INTERFACE

   !> A type useful for saving pointers to weighing functions.
   TYPE :: FmsWeightProcCfg_t
      INTEGER :: pow_value
      procedure (weight_the_field), pointer, nopass::fwf_0d_ptr   !! A pointer to the field weighting function.
      procedure (weight_the_field_1d), pointer, nopass::fwf_1d_ptr !! A pointer to the 3D field weighting function.
      procedure (weight_the_field_3d), pointer, nopass ::fwf_3d_ptr !! A pointer to the 3D field weighting function.
   CONTAINS
      procedure :: initialize => initialize_weight_proc_cfg
   END TYPE FmsWeightProcCfg_t

CONTAINS

   !> A function to linearly weight scalar fields.
   PURE REAL FUNCTION weight_the_field_0d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val        !< The field values.
      REAL, INTENT(in) :: weight           !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value     !< The weighting exponent.
      weight_the_field_0d_p1 = field_val * weight
   END FUNCTION weight_the_field_0d_p1

   !> A function to quadraticaly weight scalar fields.
   PURE REAL FUNCTION weight_the_field_0d_p2 ( field_val, weight,  pow_value  )
      REAL, INTENT(in) :: field_val        !< The field values.
      REAL, INTENT(in) :: weight           !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value    !< The weighting exponent.
      REAL :: fTw
      fTw =  field_val * weight
      weight_the_field_0d_p2 = fTw * fTw
   END FUNCTION weight_the_field_0d_p2

   !> A function to weight scalar fields by a an exponent.
   PURE REAL FUNCTION weight_the_field_0d_pp ( field_val, weight, pow_value  )
      REAL, INTENT(in) :: field_val       !< The field values.
      REAL, INTENT(in) :: weight          !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value    !< The weighting exponent.
      weight_the_field_0d_pp = (field_val * weight) ** pow_value
   END FUNCTION weight_the_field_0d_pp

   !> A function to linearly weight 1D fields.
   PURE FUNCTION weight_the_field_1d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)    !< The field values.
      REAL, INTENT(in) :: weight          !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value    !< The weighting exponent.
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d_p1
      ALLOCATE(weight_the_field_1d_p1(size(field_val)))
      DO i = 1, size(field_val)
         weight_the_field_1d_p1(i) = field_val(i) * weight
      END DO
   END FUNCTION weight_the_field_1d_p1

   !>A function to quadraticaly weight 1D fields.
   PURE FUNCTION weight_the_field_1d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)    !< The field values.
      REAL, INTENT(in) :: weight          !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value    !< The weighting exponent.
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d_p2
      !!TODO:verify the allocate
      ALLOCATE(weight_the_field_1d_p2, mold=field_val)
      DO i = 1, size(field_val)
         weight_the_field_1d_p2(i) = field_val(i) * field_val(i) * weight * weight
      END DO
   END FUNCTION weight_the_field_1d_p2

   !> A function to weight 1D fields with an exponent.
   PURE FUNCTION weight_the_field_1d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)    !< The field values.
      REAL, INTENT(in) :: weight          !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value    !< The weighting exponent.
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d_pp
      ALLOCATE(weight_the_field_1d_pp(size(field_val)))
      DO i = 1, size(field_val)
         weight_the_field_1d_pp(i) = (field_val(i) * weight) ** pow_value
      END DO
   END FUNCTION weight_the_field_1d_pp

   !> A function to linearly weight 3D fields.
   PURE FUNCTION weight_the_field_3d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:) !< The field values.
      REAL, INTENT(in) :: weight           !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value     !< The weighting exponent.
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d_p1
      ALLOCATE(weight_the_field_3d_p1, mold=field_val)
      weight_the_field_3d_p1 = field_val * weight
   END FUNCTION weight_the_field_3d_p1

   !>A function to quadraticaly weight 3D fields.
   PURE FUNCTION weight_the_field_3d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:) !< The field values.
      REAL, INTENT(in) :: weight           !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value     !< The weighting exponent.
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d_p2
      ALLOCATE(weight_the_field_3d_p2, mold=field_val)
      weight_the_field_3d_p2 = field_val * field_val * weight * weight
   END FUNCTION weight_the_field_3d_p2

   !> A function to weight 3D fields with an exponent.
   PURE FUNCTION weight_the_field_3d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)   !< The field values.
      REAL, INTENT(in) :: weight             !< The weighting coefficient.
      INTEGER, INTENT(in) :: pow_value       !< The weighting exponent.
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d_pp
      ALLOCATE(weight_the_field_3d_pp, mold=field_val)
      weight_the_field_3d_pp = (field_val * weight) ** pow_value
   END FUNCTION weight_the_field_3d_pp

   !> @brief A routine to initilize WEIGHT_PROC_CFG type.
   subroutine initialize_weight_proc_cfg ( this, pow_value )
      class (FmsWeightProcCfg_t),intent(inout) ::this !<The instance of the class that this function is bound to.
      integer, intent(in) :: pow_value  !< The weighting exponent.

      this%pow_value = pow_value

      if ( pow_value == 1) then
         this%fwf_0d_ptr => weight_the_field_0d_p1
         this%fwf_1D_ptr => weight_the_field_1d_p1
         this%fwf_3D_ptr => weight_the_field_3d_p1
      else if ( pow_value == 2 ) then
         this%fwf_0d_ptr => weight_the_field_0d_p2
         this%fwf_1d_ptr => weight_the_field_1d_p2
         this%fwf_3D_ptr => weight_the_field_3d_p2
      else
         this%fwf_0d_ptr => weight_the_field_0d_pp
         this%fwf_1D_ptr => weight_the_field_1d_pp
         this%fwf_3D_ptr => weight_the_field_3d_pp
      end if

   end subroutine initialize_weight_proc_cfg


END MODULE fms_diag_weight_procs_mod
