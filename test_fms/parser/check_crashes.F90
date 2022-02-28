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

program check_crashes
!> @brief  This programs tests if the public subroutines in parser/yaml_parser.F90
!! crash as expected
#ifdef use_yaml
use yaml_parser_mod
use mpp_mod
use fms_mod, only : fms_init, fms_end

implicit none

integer :: io_status                                  !< io_status when reading a namelist
logical :: missing_file                    = .false.  !< try to open files that do not exist
logical :: bad_conversion                  = .false.  !< try type conversions that are not possible
logical :: missing_key                     = .false.  !< try to get the value of a key that does not exist
logical :: wrong_buffer_size_key_id        = .false.  !< try to send an array of key_id that is the wrong size
logical :: wrong_buffer_size_block_id      = .false.  !< try to send an array of block_id that is the wrong size
logical :: get_key_name_bad_key_id         = .false.  !< try to send a bad key_id to get_key_name
logical :: get_key_value_bad_key_id        = .false.  !< try to send a bad key_id to get_key_value
logical :: get_block_ids_bad_id            = .false.  !< try to send a bad file_id to get_block_ids
logical :: get_key_name_bad_id             = .false.  !< try to send a bad file_id to get_key_name
logical :: get_key_value_bad_id            = .false.  !< try to send a bad file_id to get_key_value
logical :: get_num_blocks_bad_id           = .false.  !< try to send a bad file_id to get_num_blocks
logical :: get_value_from_key_bad_id       = .false.  !< try to send a bad file_id to get_value_from_key
logical :: get_nkeys_bad_id                = .false.  !< try to send a bad file_id to get_nkeys
logical :: get_key_ids_bad_id              = .false.  !< try to send a bad file_id to get_key_ids
logical :: get_key_ids_bad_block_id        = .false.  !< try to send a bad block_id to get_key_ids
logical :: get_nkeys_bad_block_id          = .false.  !< try to send a bad block_id to get_nkeys
logical :: get_block_ids_bad_block_id      = .false.  !< try to send a bad block_id to get_block_ids
logical :: get_num_blocks_bad_block_id     = .false.  !< try to send a bad block_id to get_num_blocks
logical :: get_value_from_key_bad_block_id = .false.  !< try to send a bad block_id to get_value_from_key

namelist / check_crashes_nml / missing_file, bad_conversion, missing_key, get_block_ids_bad_id, &
                               get_key_name_bad_id, get_key_value_bad_id, get_num_blocks_bad_id, &
                               get_value_from_key_bad_id, get_nkeys_bad_id, get_key_ids_bad_id, &
                               get_key_name_bad_key_id, get_key_value_bad_key_id, get_key_ids_bad_block_id, &
                               get_nkeys_bad_block_id, get_block_ids_bad_block_id,  get_num_blocks_bad_block_id, &
                               get_value_from_key_bad_block_id, &
                               wrong_buffer_size_key_id, wrong_buffer_size_block_id

call fms_init

read (input_nml_file, check_crashes_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>check_crashes: Error reading input.nml')

if (missing_file)                    call check_read_and_parse_file_missing
if (get_block_ids_bad_id)            call check_get_block_ids_bad_id
if (get_key_name_bad_id)             call check_get_key_name_bad_id
if (get_key_value_bad_id)            call check_get_key_value_bad_id
if (get_num_blocks_bad_id)           call check_get_num_blocks_bad_id
if (get_value_from_key_bad_id)       call check_get_value_from_key_bad_id
if (get_nkeys_bad_id)                call check_get_nkeys_bad_id
if (get_key_ids_bad_id)              call check_get_key_ids_bad_id
if (bad_conversion)                  call check_bad_conversion
if (missing_key)                     call check_missing_key
if (wrong_buffer_size_key_id)        call check_wrong_buffer_size_key_id
if (wrong_buffer_size_block_id)      call check_wrong_buffer_size_block_id
if (get_key_name_bad_key_id)         call check_get_key_name_bad_key_id
if (get_key_value_bad_key_id)        call check_get_key_value_bad_key_id
if (get_key_ids_bad_block_id)        call check_get_key_ids_bad_block_id
if (get_nkeys_bad_block_id)          call check_get_nkeys_bad_block_id
if (get_block_ids_bad_block_id)      call check_get_block_ids_bad_block_id
if (get_num_blocks_bad_block_id)     call check_get_num_blocks_bad_block_id
if (get_value_from_key_bad_block_id) call check_get_value_from_key_bad_block_id

call fms_end

contains
!> @brief This is to check if the parser crashes correctly if user sends a bad block_id to get_key_ids
subroutine check_get_key_ids_bad_block_id
   integer :: yaml_file_id !< file_id for a yaml file
   integer :: key_ids(10)  !< array of key ids

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_key_ids (yaml_file_id, -40, key_ids)

end subroutine check_get_key_ids_bad_block_id

!> @brief This is to check if the parser crashes correctly if user sends a bad block_id to get_nkeys
subroutine check_get_nkeys_bad_block_id
   integer :: yaml_file_id !< file_id for a yaml file
   integer :: nkeys        !< number of keys

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   nkeys = get_nkeys(yaml_file_id, 9999)

end subroutine check_get_nkeys_bad_block_id

!> @brief This is to check if the parser crashes correctly if user sends a bad parent_block_id to get_block_ids
subroutine check_get_block_ids_bad_block_id
   integer :: yaml_file_id !< file_id for a yaml file
   integer :: block_ids(10)!< array of block ids

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_block_ids(yaml_file_id, "varList", block_ids, parent_block_id=-40)

end subroutine check_get_block_ids_bad_block_id

!> @brief This is to check if the parser crashes correctly if user sends a bad parent_block_id to get_num_blocks
subroutine check_get_num_blocks_bad_block_id
   integer :: yaml_file_id !< file_id for a yaml file
   integer :: nblocks      !< number of blocks

   yaml_file_id = open_and_parse_file("diag_table.yaml")

   nblocks = get_num_blocks(yaml_file_id, "varList", parent_block_id=-30)

end subroutine check_get_num_blocks_bad_block_id

!> @brief This is to check if the parser crashes correctly if user sends a bad parent_block_id to get_value_from_key
subroutine check_get_value_from_key_bad_block_id
   integer :: yaml_file_id !< file_id for a yaml file
   integer :: key_value    !< integer buffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_value_from_key(yaml_file_id, 999, "mullions", key_value)

end subroutine check_get_value_from_key_bad_block_id

!> @brief This is to check if the parser crashes correctly if user tries to open a missing file.
subroutine check_read_and_parse_file_missing
   integer :: yaml_file_id !< file_id for a yaml file

   yaml_file_id = open_and_parse_file("missing")
end subroutine check_read_and_parse_file_missing

!> @brief This is to check if the parser crashes correctly if user sends an invalid file id to get_block_ids
subroutine check_get_block_ids_bad_id
   integer :: block_ids(10) !< array of block ids

   call get_block_ids(-40, "diagFiles", block_ids)
end subroutine check_get_block_ids_bad_id

!> @brief This is to check if the parser crashes correctly if user sends an invalid file id to get_key_name
subroutine check_get_key_name_bad_id
   character(len=10) :: buffer !< string buffer

   call get_key_name(-45, 1, buffer)
end subroutine check_get_key_name_bad_id

!> @brief This is to check if the parser crashes correctly if user sends an invalid file id to get_key_value
subroutine check_get_key_value_bad_id
   character(len=10) :: buffer !< string buffer

   call get_key_value(-45, 1, buffer)
end subroutine check_get_key_value_bad_id

!> @brief This is to check if the parser crashes correctly if user sends an invalid file id to get_num_blocks
subroutine check_get_num_blocks_bad_id
   integer :: nblocks !< number of blocks

   nblocks = get_num_blocks(-45, "diagFiles")
end subroutine check_get_num_blocks_bad_id

!> @brief This is to check if the parser crashes correctly if user sends an invalid file id to get_value_from_key
subroutine check_get_value_from_key_bad_id
   character(len=10) :: string_buffer !< string buffer

   call get_value_from_key(-45, 1, "varName", string_buffer)
end subroutine check_get_value_from_key_bad_id

!> @brief This is to check if the parser crashes correctly if user sends an invalid file id to get_nkeys
subroutine check_get_nkeys_bad_id
   integer :: nkeys !< number of keys

   nkeys = get_nkeys(-45, 1)
end subroutine check_get_nkeys_bad_id

!> @brief This is to check if the parser crashes correctly if user sends an invalid file id to get_key_ids
subroutine check_get_key_ids_bad_id
   integer :: key_ids(10) !< array of key ids

   call get_key_ids(-45, 1, key_ids)
end subroutine check_get_key_ids_bad_id

!> @brief This is to check if the parser crashes correctly if user sends a buffer of the wrong type
subroutine check_bad_conversion
   integer :: yaml_file_id !< file_id for a yaml file
   real :: buffer          !< real buffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_value_from_key(yaml_file_id, 9, "varName", buffer)
end subroutine check_bad_conversion

!> @brief This is to check if the parser crashes correctly if user tries to get they value for a key
!! that doesn't exist
subroutine check_missing_key
   integer :: yaml_file_id !< file_id for a yaml file
   real    :: buffer       !< string bufffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_value_from_key(yaml_file_id, 9, "missing", buffer)
end subroutine check_missing_key

!> @brief This is to check if the parser crashes correctly if user sends an invalid key id to get_key_name
subroutine check_get_key_name_bad_key_id
    integer :: yaml_file_id     !< file_id for a yaml file
    character(len=10) :: buffer !< string buffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_key_name(yaml_file_id, 666, buffer)

end subroutine check_get_key_name_bad_key_id

!> @brief This is to check if the parser crashes correctly if user sends an invalid key id to get_key_value
subroutine check_get_key_value_bad_key_id
    integer :: yaml_file_id     !< file_id for a yaml file
    character(len=10) :: buffer !< string buffer

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_key_value(yaml_file_id, 666, buffer)

end subroutine check_get_key_value_bad_key_id

!> @brief This is to check if the parser crashes correctly if user sends an a key_id array that is that the correct
!! size to get_key_ids
subroutine check_wrong_buffer_size_key_id
   integer :: yaml_file_id !< file_id for a yaml file
   integer :: key_ids(1)   !< array of key ids

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_key_ids(yaml_file_id, 19, key_ids)

end subroutine check_wrong_buffer_size_key_id

!> @brief This is to check if the parser crashes correctly if user sends an a block_id array that is that the correct
!! size to get_block_ids
subroutine check_wrong_buffer_size_block_id
   integer :: yaml_file_id !< file_id for a yaml file
   integer :: block_ids(10)!< array of block ids

   yaml_file_id = open_and_parse_file("diag_table.yaml")
   call get_block_ids(yaml_file_id, "diag_files", block_ids)

end subroutine check_wrong_buffer_size_block_id
#endif
end program check_crashes
