module fms_diag_yaml_mod

use diag_data_mod, only: diag_files_type, diag_fields_type

integer, parameter :: basedate_size = 7

!> Object that holds the information of the diag_yaml
type diag_yaml_object
 character(len=:), allocatable, private :: diag_title                   !< Experiment name
 integer, private, dimension (basedate_size) :: diag_basedate           !< basedate array
 type(diag_files_type), allocatable, private, dimension (:) :: diag_files!< History file info
 type(diag_fields_type), allocatable, private, dimension (:,:) :: diag_fields !< Diag fields info
 contains
 procedure :: title => get_title        !< Returns the title
 procedure :: basedate => get_basedate  !< Returns the basedate array
end type diag_yaml_object
type (diag_yaml_object) :: diag_yaml

public :: get_title, get_basedate

contains

!> \brief Returns the basedate as an integer array
pure function get_basedate (diag_yaml) result (diag_basedate)
class (diag_yaml_object), intent(in) :: diag_yaml               !< The diag_yaml
integer, dimension (basedate_size) :: diag_basedate !< Basedate array result to return
diag_basedate = diag_yaml%diag_basedate
end function get_basedate
!> \brief Returns the title of the diag table as an allocated string
pure function get_title (diag_yaml) result (diag_title)
class (diag_yaml_object), intent(in) :: diag_yaml      !< The diag_yaml
character(len=:),allocatable :: diag_title !< Basedate array result to return
 diag_title = diag_yaml%diag_title
end function get_title

!> \brief Compares two field type variables
pure logical function is_field_type_null (in1)
type(diag_fields_type), intent(in) :: in1
is_field_type_null = .true.
end function is_field_type_null

!!TODO
!> \brief looks for a diag_field based on it's name.
!! Returns null if field is not found.
!type(diag_fields_type)function get_diag_table_field (field_name) result (field)
! character(len=*), intent(IN) :: field_name
! integer :: i
! do i = 1,size(diag_fields)
!     if (trim(field_name) == trim(fms_c2f_string(diag_fields(i)%fname))) then
!          field = diag_fields(i)
!write (6,*) field_name//" Found"
!
!          return
!     endif
! enddo
! field = null_field_type
!
!end function get_diag_table_field



end module fms_diag_yaml_mod
