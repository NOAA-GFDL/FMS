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

!> @brief  This programs tests the public subroutines in parser/yaml_parser.F90
program test_read_and_parse_file

#ifdef use_yaml
use yaml_parser_mod
use mpp_mod
use fms_mod, only : fms_init, fms_end
use platform_mod

implicit none

integer               :: yaml_file_id1   !< file id of a yaml file
integer               :: nfiles          !< number of files
integer               :: nvariables      !< number of variables
integer, allocatable  :: file_ids(:)     !< array of file ids
integer, allocatable  :: variable_ids(:) !< array of variable ids
integer               :: yaml_file_id2   !< file id of a yaml file
integer               :: nentries        !< number of entries
integer, allocatable  :: entries_ids(:)  !< array of entries ids
integer               :: zero            !< dummy integer buffer
character(len=20)     :: string_buffer   !< string buffer
integer(kind=i4_kind) :: i4_buffer       !< i4 buffer
integer(kind=i8_kind) :: i8_buffer       !< i8 buffer
real(kind=r4_kind)    :: r4_buffer       !< r4 buffer
real(kind=r8_kind)    :: r8_buffer       !< r8 buffer
integer               :: nkeys           !< number of keys
integer, allocatable  :: key_ids(:)      !< array of key ids
character(len=20)     :: key_name        !< the name of the key
character(len=20)     :: key_value       !< the value of a key
logical               :: logical_buffer  !< logical buffer

call fms_init

!< Test open_and_parse_file
yaml_file_id1 = open_and_parse_file("diag_table.yaml")
if (yaml_file_id1 .ne. 0) call mpp_error(FATAL, "The yaml_file_id for this file should be 0")

!< Test if multiple files can be opened
yaml_file_id2 = open_and_parse_file("data_table.yaml")
if (yaml_file_id2 .ne. 1) call mpp_error(FATAL, "The yaml_file_id for this file should be 1")

!< -----------------------------------

!< Test get_num_blocks
nfiles = get_num_blocks(yaml_file_id1, "diag_files")
if (nfiles .ne. 2) call mpp_error(FATAL, "There should be only 2 diag_files")

!< Test if a different yaml file id will work
nentries = get_num_blocks(yaml_file_id2, "data_table")
if (nentries .ne. 2) call mpp_error(FATAL, "There should be only 2 entries")

!< Try to look for a block that does not exist!
zero = get_num_blocks(yaml_file_id2, "diag_files")
if (zero .ne. 0) call mpp_error(FATAL, "'diag_files' should not exist in this file")

!< Try the parent block_id optional argument
nvariables = get_num_blocks(yaml_file_id1, "varlist", parent_block_id=3) !< Number of variables that belong to the atmos_daily file in the diag_table.yaml
if (nvariables .ne. 2) call mpp_error(FATAL, "There should only be 2 variables in the atmos_daily file")

!< -----------------------------------

!< Test get_block_ids
allocate(file_ids(nfiles))
call get_block_ids(yaml_file_id1, "diag_files", file_ids)
if(file_ids(1) .ne. 3 .or. file_ids(2) .ne. 21) call mpp_error(FATAL, "The file_ids are wrong!")

!< Test to see if a diffrent yaml file id will work
allocate(entries_ids(nentries))
call get_block_ids(yaml_file_id2, "data_table", entries_ids)
if(entries_ids(1) .ne. 1 .or. entries_ids(2) .ne. 8) call mpp_error(FATAL, "The entry_ids are wrong!")

!< Try the parent block id optional argument
allocate(variable_ids(nvariables))
call get_block_ids(yaml_file_id1, "varlist", variable_ids, parent_block_id=3)
if (variable_ids(1) .ne. 9 .or. variable_ids(2) .ne. 15) call mpp_error(FATAL, "The variable_ids are wrong!")

!< Error check: *_ids is not the correct size

!< -----------------------------------

!< Test get_value_from_key
!! Try get_value_from_key using a string buffer
call get_value_from_key(yaml_file_id1, variable_ids(1), "varName", string_buffer)
if (trim(string_buffer) .ne. "tdata") call mpp_error(FATAL, "varName was not read correctly!")

!! Try get_value_from_key using a i4 buffer
call get_value_from_key(yaml_file_id1, variable_ids(1), "mullions", i4_buffer)
if (i4_buffer .ne. int(10, kind=i4_kind)) call mpp_error(FATAL, "mullions was not read correctly as an i4!")

!! Try get_value_from_key using a i8 buffer
call get_value_from_key(yaml_file_id1, variable_ids(1), "mullions", i8_buffer)
if (i8_buffer .ne. int(10, kind=i8_kind)) call mpp_error(FATAL, "mullions was not read correctly as an i8!")

!! Try get_value_from_key using a r4 buffer
call get_value_from_key(yaml_file_id1, variable_ids(1), "fill_value", r4_buffer)
if (r4_buffer .ne. real(-999.9, kind=r4_kind)) call mpp_error(FATAL, "fill_value was not read correctly as an r4!")

!! Try get_value_from_key using a r8 buffer
call get_value_from_key(yaml_file_id1, variable_ids(1), "fill_value", r8_buffer)
if (abs(r8_buffer - real(-999.9, kind=r8_kind)) .gt. 5d-5) call mpp_error(FATAL, "fill_value was not read correctly as an r8!")

!! Try get_value_from_key using a logical buffer
logical_buffer = .true.
call get_value_from_key(yaml_file_id2, entries_ids(2), "do_data_bug", logical_buffer)
if (logical_buffer) call mpp_error(FATAL, "do_data_bug was not read correctly as a logical")

logical_buffer = .false.
call get_value_from_key(yaml_file_id2, entries_ids(2), "use_data_bug", logical_buffer)
if (.not. logical_buffer) call mpp_error(FATAL, "use_data_bug was not read correctly as a logical")

!! Try the is_optional argument on an key that does not exist
string_buffer = ""
call get_value_from_key(yaml_file_id1, variable_ids(1), "NANANANA", string_buffer, is_optional=.true.)
if (trim(string_buffer) .ne. "") call mpp_error(FATAL, "string_buffer was set when they key does not exist?")

!< -----------------------------------

!< Test nkeys
nkeys = get_nkeys(yaml_file_id1, variable_ids(1))
if (nkeys .ne. 5) call mpp_error(FATAL, "The number of keys was not read correctly")

!< -----------------------------------

!< Test get_key_ids
allocate(key_ids(nkeys))
call get_key_ids(yaml_file_id1, variable_ids(1), key_ids)
if (key_ids(1) .ne. 10 .or. key_ids(2) .ne. 11 .or. key_ids(3) .ne. 12 .or. key_ids(4) .ne. 13 .or. key_ids(5) .ne. 14) call mpp_error(FATAL, "The key ids obtained are wrong")

!< -----------------------------------

!< Test get_key_name
call get_key_name(yaml_file_id1, key_ids(1), key_name)
if ((trim(key_name) .ne. "varName")) call mpp_error(FATAL, "get_key_name did not output the correct name")

!< Test get_key_value
call get_key_value(yaml_file_id1, key_ids(1), key_value)
if ((trim(key_value) .ne. "tdata")) call mpp_error(FATAL, "get_key_name did not output the correct name")

deallocate(key_ids)
deallocate(variable_ids)
deallocate(entries_ids)
deallocate(file_ids)

call fms_end
#endif
end program
