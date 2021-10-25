program check_crashes
#ifdef use_yaml
use yaml_parser_mod
use mpp_mod
use fms_mod, only : fms_init, fms_end

implicit none

integer :: io_status
logical :: missing_file        = .false.
logical :: bad_conversion      = .false.
logical :: missing_key         = .false.
logical :: wrong_buffer_size   = .false.
logical :: get_key_name_bad_id = .false.

namelist / check_crashes_nml / missing_file, bad_conversion, missing_key, wrong_buffer_size, get_key_name_bad_id

call fms_init

read (input_nml_file, check_crashes_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>check_crashes: Error reading input.nml')

!< Bad file id
if (missing_file)        call check_read_and_parse_file_missing
if (bad_conversion)      call check_bad_conversion
if (missing_key)         call check_missing_key
if (wrong_buffer_size)   call check_wrong_buffer_size
if (get_key_name_bad_id) call check_get_key_name_bad_id

call fms_end

contains

!> @brief This is to check if the parser crashes correctly if user tries to open a missing file. 
subroutine check_read_and_parse_file_missing
   integer :: yaml_file_id
   yaml_file_id = open_and_parse_file("missing")
end subroutine check_read_and_parse_file_missing

!> @brief This is to check if the parser crashes correctly if user sends a buffer of the wrong type
subroutine check_bad_conversion
   integer :: yaml_file_id
   real :: buffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_value_from_key(yaml_file_id, 9, "varName", buffer)
end subroutine check_bad_conversion

!> @brief This is to check if the parser crashes correctly if user tries to get they value for a key
!! that doesn't exist
subroutine check_missing_key
   integer :: yaml_file_id
   real :: buffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_value_from_key(yaml_file_id, 9, "missing", buffer)
end subroutine check_missing_key

subroutine check_get_key_name_bad_id
    integer :: yaml_file_id
    character(len=10) :: buffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_key_value(yaml_file_id, 666, buffer)

end subroutine check_get_key_name_bad_id

subroutine check_wrong_buffer_size
   integer :: yaml_file_id
   integer :: file_ids(1)

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_block_ids(yaml_file_id, "diag_files", file_ids)
   print *, file_ids

end subroutine check_wrong_buffer_size
#endif
end program check_crashes
