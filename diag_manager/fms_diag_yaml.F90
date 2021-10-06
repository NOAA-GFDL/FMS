module fms_diag_yaml_mod

use fms_diag_data_mod, only: diag_files_type, diag_fields_type

contains
!> \brief Compares two field type variables
pure logical function is_field_type_null (in1)
type(diag_fields_type), intent(in) :: in1
is_field_type_null = (in1%ikind == DIAG_NULL)
end function is_field_type_null

!!TODO
!> \brief looks for a diag_field based on it's name.
!! Returns null if field is not found.
type(diag_fields_type)function get_diag_table_field (field_name) result (field)
 character(len=*), intent(IN) :: field_name
 integer :: i
! do i = 1,size(diag_fields)
!     if (trim(field_name) == trim(fms_c2f_string(diag_fields(i)%fname))) then
!          field = diag_fields(i)
!write (6,*) field_name//" Found"
!
!          return
!     endif
! enddo
! field = null_field_type

end function get_diag_table_field



end module fms_diag_yaml_mod
