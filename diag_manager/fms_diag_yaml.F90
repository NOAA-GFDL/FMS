module fms_diag_yaml_mod

use fms_diag_yaml_object_mod, only: diag_yaml_files_type, diag_yaml_files_var_type, diag_yaml_files_obj_init
use yaml_parser_mod

implicit none

integer, parameter :: basedate_size = 6

!> Object that holds the information of the diag_yaml
type diag_yaml_object
 character(len=:), allocatable, private :: diag_title                   !< Experiment name
 integer, private, dimension (basedate_size) :: diag_basedate           !< basedate array
 type(diag_yaml_files_type), allocatable, private, dimension (:) :: diag_files!< History file info
 type(diag_yaml_files_var_type), allocatable, private, dimension (:) :: diag_fields !< Diag fields info
end type diag_yaml_object

type (diag_yaml_object) :: diag_yaml

contains

subroutine diag_yaml_object_init
integer              :: diag_yaml_id   !< Id for the diag_table yaml
integer              :: nfiles         !< Number of files in the diag_table yaml
integer, allocatable :: file_ids(:)    !< Ids of the files in the diag_table yaml

integer :: i
integer :: total_nvars
type(diag_yaml_files_type) :: temp_file

diag_yaml_id = open_and_parse_file("diag_table.yaml")

call get_value_from_key(diag_yaml_id, 0, "title", diag_yaml%diag_title)
call get_value_from_key(diag_yaml_id, 0, "baseDate", diag_yaml%diag_basedate)

nfiles = get_num_blocks(diag_yaml_id, "diag_files")
allocate(diag_yaml%diag_files(nfiles))
allocate(file_ids(nfiles))
call get_block_ids(diag_yaml_id, "diag_files", file_ids)

total_nvars = get_total_num_vars(diag_yaml_id, file_ids)
allocate(diag_yaml%diag_fields(total_nvars))

do i = 1, nfiles
   call init_diag_yaml_files_obj_init(temp_file)
   call get_value_from_key(diag_yaml_id, file_ids(i), "file_name", temp_file%file_fname)
   call get_value_from_key(diag_yaml_id, file_ids(i), "freq_units", temp_file%file_frequnit)
   call get_value_from_key(diag_yaml_id, file_ids(i), "freq", temp_file%file_freq)
   call get_value_from_key(diag_yaml_id, file_ids(i), "unlimdim", temp_file%file_timeunit)
   call get_value_from_key(diag_yaml_id, file_ids(i), "time_units", temp_file%file_unlimdim)
   call get_value_from_key(diag_yaml_id, file_ids(i), "write_file", temp_file%string_file_write, is_optional=.true.)
enddo

deallocate(file_ids)
end subroutine

function get_total_num_vars(diag_yaml_id, file_ids) &
result(total_nvars)

integer, intent(in) :: diag_yaml_id   !< Id for the diag_table yaml
integer, intent(in) :: file_ids(:)    !< Ids of the files in the diag_table yaml
integer :: total_nvars
integer :: buffer

integer :: i

total_nvars = 0
do i = 1, size(file_ids,1)
    buffer = get_num_blocks(diag_yaml_id, "varlist", parent_block_id=file_ids(i))
    total_nvars = total_nvars + buffer
end do

end function

end module fms_diag_yaml_mod
