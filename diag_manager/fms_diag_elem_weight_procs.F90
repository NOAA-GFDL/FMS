MODULE fms_diag_elem_weight_procs_mod
   USE platform_mod

   implicit none
   INTERFACE addwf
      module procedure addwf_r4
      module procedure addwf_r8
      module procedure addwf_i4
      module procedure addwf_i8
   END INTERFACE

CONTAINS

   ELEMENTAL PURE SUBROUTINE addwf_r4(buff,  field, weight, pow_value )
      REAL(r4_kind), INTENT(inout) :: buff
      REAL(r4_kind), INTENT(IN) :: field
      REAL, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value

      SELECT  CASE(pow_value)
       CASE (1)
         buff = buff + weight * field
       CASE (2)
         buff = buff + (weight * field) *  (weight * field)
       CASE  default
         buff = buff + (weight * field) ** pow_value
      END SELECT
   END SUBROUTINE addwf_r4

   ELEMENTAL PURE SUBROUTINE addwf_r8(buff,  field, weight, pow_value )
      REAL(r8_kind), INTENT(inout) :: buff
      REAL(r8_kind) ,INTENT(IN) :: field
      REAL, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value

      SELECT  CASE(pow_value)
       CASE (1)
         buff = buff + weight * field
       CASE (2)
         buff = buff + (weight * field) *  (weight * field)
       CASE  default
         buff = buff + (weight * field) ** pow_value
      END SELECT
   END SUBROUTINE addwf_r8

   ELEMENTAL PURE SUBROUTINE addwf_i4(buff,  field, weight, pow_value )
      INTEGER(i4_kind), INTENT(inout) :: buff
      INTEGER(i4_kind), INTENT(IN) :: field
      INTEGER, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value
      SELECT  CASE(pow_value)
       CASE (1)
         buff = buff + weight * field
       CASE (2)
         buff = buff + (weight * field) *  (weight * field)
       CASE  default
         buff = buff + (weight * field) ** pow_value
      END SELECT
   END SUBROUTINE addwf_i4

   ELEMENTAL PURE SUBROUTINE addwf_i8(buff,  field, weight, pow_value )
      INTEGER(i8_kind), INTENT(inout) :: buff
      INTEGER(i8_kind) ,INTENT(IN) :: field
      INTEGER, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value

      SELECT  CASE(pow_value)
       CASE (1)
         buff = buff + weight * field
       CASE (2)
         buff = buff + (weight * field) *  (weight * field)
       CASE  default
         buff = buff + (weight * field) ** pow_value
      END SELECT
   END SUBROUTINE addwf_i8
END MODULE fms_diag_elem_weight_procs_mod

