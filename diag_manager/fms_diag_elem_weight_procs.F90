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

   ELEMENTAL PURE  REAL(r4_kind) FUNCTION addwf_r4(buff,  field, weight, pow_value )
      REAL(r4_kind), INTENT(in) :: buff
      REAL(r4_kind), INTENT(IN) :: field
      REAL, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value

      SELECT  CASE(pow_value)
       CASE (1)
        addwf_r4 = buff + weight * field
       CASE (2)
        addwf_r4 = buff + (weight * field) *  (weight * field)
       CASE  default
        addwf_r4 = buff + (weight * field) ** pow_value
      END SELECT
   END FUNCTION addwf_r4

   ELEMENTAL PURE  REAL(r8_kind) FUNCTION addwf_r8(buff,  field, weight, pow_value )
      REAL(r8_kind), INTENT(in) :: buff
      REAL(r8_kind) ,INTENT(IN) :: field
      REAL, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value

      SELECT  CASE(pow_value)
       CASE (1)
        addwf_r8 = buff + weight * field
       CASE (2)
        addwf_r8 = buff + (weight * field) *  (weight * field)
       CASE  default
        addwf_r8 = buff + (weight * field) ** pow_value
      END SELECT
   END FUNCTION addwf_r8

   ELEMENTAL PURE INTEGER(i4_kind) FUNCTION addwf_i4(buff,  field, weight, pow_value )
      INTEGER(i4_kind), INTENT(in) :: buff
      INTEGER(i4_kind), INTENT(IN) :: field
      INTEGER, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value
      SELECT  CASE(pow_value)
       CASE (1)
        addwf_i4 = buff + weight * field
       CASE (2)
        addwf_i4 = buff + (weight * field) *  (weight * field)
       CASE  default
        addwf_i4 = buff + (weight * field) ** pow_value
      END SELECT
   END FUNCTION addwf_i4

   ELEMENTAL PURE INTEGER(i8_kind) FUNCTION addwf_i8(buff,  field, weight, pow_value )
      INTEGER(i8_kind), INTENT(in) :: buff
      INTEGER(i8_kind) ,INTENT(IN) :: field
      INTEGER, INTENT(IN) ::  weight
      INTEGER, INTENT(IN) :: pow_value

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

