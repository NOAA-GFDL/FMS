!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @defgroup fms_diag_elem_weight_procs_mod fms_diag_elem_weight_procs_mod
!> @ingroup diag_manager
!> @brief fms_diag_elem_weight_procs_mod Contains elemental functions for uddating
!! one element of a buffer array with field data.
!!
!> @author Miguel Zuniga
!!
!! <TT>fms_diag_elem_weight_procs_mod</TT> Contains elemental functions for uddating
!! one element of a buffer array with field data,
!!
!> @file
!> @brief File for @ref fms_diag_elem_weight_procs_mod
!> @addtogroup fms_diag_elem_weight_procs_mod
!> @{
MODULE fms_diag_elem_weight_procs_mod
   USE platform_mod

   implicit none

  !> @brief Interface for the elemental function addwf, which
  !! Calculates and returns the value given by this formula:
  !! returned_value  = buff  + (weight * field)**pow_value
  !> @ingroup fms_diag_elem_weight_procs_mod
   INTERFACE addwf
      module procedure addwf_r4
      module procedure addwf_r8
      module procedure addwf_i4
      module procedure addwf_i8
   END INTERFACE

CONTAINS

   !!TODO: Note that in the functions below, the case for pow_value == 2 was
   !! not in the original send_data_3d code and the power function was used.
   !! So this case may need to be deleted if reproducability is an issue.

   !!TODO: (MDM) Discuss whether or not the pow_value should be allowed to
   !! also be real though legacy interface has it satic.

  !> @brief Calculates and returns the value given by this formula:
  !! returned_value  = buff  + (weight * field)**pow_value
  !! Special cases when pow_value is equal to 1 or 2 do not explicitly use the power function.
   ELEMENTAL REAL(r4_kind) FUNCTION addwf_r4(buff,  field, weight, pow_value )
      REAL(r4_kind), INTENT(in) :: buff !< The buffer cell (point) value
      REAL(r4_kind), INTENT(IN) :: field !< The field value
      REAL(r4_kind), INTENT(IN) ::  weight !< The weight factor for the field
      INTEGER, INTENT(IN) :: pow_value !< The power value for the power function

      SELECT  CASE(pow_value)
       CASE (1)
        addwf_r4 = buff + weight * field
       CASE (2)
        addwf_r4 = buff + (weight * field) *  (weight * field)
       CASE  default
        addwf_r4 = buff + (weight * field) ** pow_value
      END SELECT
   END FUNCTION addwf_r4

  !> @brief Calculates and returns the value given by this formula:
  !! returned_value  = buff  + (weight * field)**pow_value
  !! Special cases when pow_value is equal to 1 or 2 do not explicitly use the power function.
   ELEMENTAL REAL(r8_kind) FUNCTION addwf_r8(buff,  field, weight, pow_value )
      REAL(r8_kind), INTENT(in) :: buff !< The buffer cell (point) value
      REAL(r8_kind) ,INTENT(IN) :: field !< The field value
      REAL(r8_kind), INTENT(IN) :: weight !< The weight factor for the field
      INTEGER, INTENT(IN) :: pow_value !< The power value for the power function

      SELECT  CASE(pow_value)
       CASE (1)
        addwf_r8 = buff + weight * field
       CASE (2)
        addwf_r8 = buff + (weight * field) *  (weight * field)
       CASE  default
        addwf_r8 = buff + (weight * field) ** pow_value
      END SELECT
   END FUNCTION addwf_r8

  !> @brief Calculates and returns the value given by this formula:
  !! returned_value  = buff  + (weight * field)**pow_value
  !! Special cases when pow_value is equal to 1 or 2 do not explicitly use the power function.
   ELEMENTAL INTEGER(i4_kind) FUNCTION addwf_i4(buff,  field, weight, pow_value )
      INTEGER(i4_kind), INTENT(in) :: buff !< The buffer cell (point) value
      INTEGER(i4_kind), INTENT(IN) :: field !< The field value
      INTEGER, INTENT(IN) ::  weight !< The weight factor for the field
      INTEGER, INTENT(IN) :: pow_value !< The power value for the power function
      SELECT  CASE(pow_value)
       CASE (1)
        addwf_i4 = buff + weight * field
       CASE (2)
        addwf_i4 = buff + (weight * field) *  (weight * field)
       CASE  default
        addwf_i4 = buff + (weight * field) ** pow_value
      END SELECT
   END FUNCTION addwf_i4

  !> @brief Calculates and returns the value given by this formula:
  !! returned_value  = buff  + (weight * field)**pow_value
  !! Special cases when pow_value is equal to 1 or 2 do not explicitly use the power function.
   ELEMENTAL INTEGER(i8_kind) FUNCTION addwf_i8(buff,  field, weight, pow_value )
      INTEGER(i8_kind), INTENT(in) :: buff !< The buffer cell (point) value
      INTEGER(i8_kind) ,INTENT(IN) :: field !< The field value
      INTEGER, INTENT(IN) ::  weight !< The weight factor for the field
      INTEGER, INTENT(IN) :: pow_value !< The power value for the power function

      SELECT  CASE(pow_value)
       CASE (1)
        addwf_i8 = buff + weight * field
       CASE (2)
        addwf_i8 = buff + (weight * field) *  (weight * field)
       CASE  default
        addwf_i8 = buff + (weight * field) ** pow_value
      END SELECT
   END FUNCTION addwf_i8
END MODULE fms_diag_elem_weight_procs_mod
!> @}
! close documentation grouping

