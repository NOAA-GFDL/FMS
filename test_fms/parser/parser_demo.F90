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

program parser_demo
!> @brief  This programs demostrates how to use the parser

#ifdef use_yaml
use FMS_mod, only: fms_init, fms_end
use yaml_parser_mod
use platform_mod

implicit none

integer              :: diag_yaml_id   !< Id for the diag_table yaml
integer              :: nfiles         !< Number of files in the diag_table yaml
integer, allocatable :: file_ids(:)    !< Ids of the files in the diag_table yaml
integer              :: nvariables     !< Number of variables in the diag_table yaml
integer, allocatable :: var_ids(:)     !< Ids of the variables in the diag_table yaml
integer              :: i, j           !< For do loops
character(len=255)   :: string_buffer  !< Buffer to read strings to
integer              :: int_buffer     !< Buffer to read integers to
real(kind=r8_kind)   :: r8_buffer      !< Buffer to read r8 to

call fms_init

diag_yaml_id = open_and_parse_file("diag_table.yaml")
print *, ""

call get_value_from_key(diag_yaml_id, 0, "title", string_buffer)
print *, "title:", trim(string_buffer)

call get_value_from_key(diag_yaml_id, 0, "baseDate", string_buffer)
print *, "baseDate:", trim(string_buffer)

nfiles = get_num_blocks(diag_yaml_id, "diag_files")
allocate(file_ids(nfiles))
call get_block_ids(diag_yaml_id, "diag_files", file_ids)
print *, ""

do i = 1, nfiles
   print *, "File number:", i

   call get_value_from_key(diag_yaml_id, file_ids(i), "fileName", string_buffer)
   print *, "fileName:", trim(string_buffer)

   call get_value_from_key(diag_yaml_id, file_ids(i), "freq", int_buffer)
   print *, "freq:", int_buffer

   call get_value_from_key(diag_yaml_id, file_ids(i), "frequnit", string_buffer)
   print *, "frequnit:", trim(string_buffer)

   call get_value_from_key(diag_yaml_id, file_ids(i), "timeunit", string_buffer)
   print *, "timeunit:", trim(string_buffer)

   call get_value_from_key(diag_yaml_id, file_ids(i), "unlimdim", string_buffer)
   print *, "unlimdim:", trim(string_buffer)

   !< The number of variables that are part of the current file
   nvariables = get_num_blocks(diag_yaml_id, "varlist", parent_block_id=file_ids(i))
   allocate(var_ids(nvariables))
   call get_block_ids(diag_yaml_id, "varlist", var_ids, parent_block_id=file_ids(i))

   do j = 1, nvariables
       print *, "  Variable number:", j

       call get_value_from_key(diag_yaml_id, var_ids(j), "varName", string_buffer)
       print *, "  varName:", trim(string_buffer)

       string_buffer = ""
       call get_value_from_key(diag_yaml_id, var_ids(j), "reduction", string_buffer)
       print *, "  reduction:", trim(string_buffer)

       string_buffer = ""
       call get_value_from_key(diag_yaml_id, var_ids(j), "module", string_buffer)
       print *, "  module:", trim(string_buffer)

       r8_buffer = 0.
       call get_value_from_key(diag_yaml_id, var_ids(j), "fill_value", r8_buffer, is_optional=.true.)
       print *, "  fill_value:", r8_buffer

       string_buffer = ""
       call get_value_from_key(diag_yaml_id, var_ids(j), "outName", string_buffer, is_optional=.true.)
       print *, "  outName:", trim(string_buffer)

       string_buffer = ""
       call get_value_from_key(diag_yaml_id, var_ids(j), "kind", string_buffer, is_optional=.true.)
       print *, "  kind:", trim(string_buffer)

       int_buffer = 0.
       call get_value_from_key(diag_yaml_id, var_ids(j), "mullions", int_buffer, is_optional=.true.)
       print *, "  mullions:", int_buffer

       print *, ""
   end do
   deallocate(var_ids)
   print *, ""
enddo
deallocate(file_ids)
call fms_end

#endif
end program parser_demo
