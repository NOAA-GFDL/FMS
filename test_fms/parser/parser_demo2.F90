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
integer              :: i, j, k        !< For do loops
integer              :: nkeys          !< Number of keys
integer, allocatable :: key_ids(:)     !< Ids of keys in the diag_table_yaml
character(len=255)   :: key_value      !< The value of a key
character(len=255)   :: key_name       !< The name of a key

call fms_init

diag_yaml_id = open_and_parse_file("diag_table.yaml")
print *, ""

nkeys = get_nkeys(diag_yaml_id, 0)
allocate(key_ids(nkeys))
call get_key_ids(diag_yaml_id, 0, key_ids)

do i = 1, nkeys
    call get_key_name(diag_yaml_id, key_ids(i), key_name)
    call get_key_value(diag_yaml_id, key_ids(i), key_value)
    print *, "Key:", trim(key_name), " Value:", trim(key_value)
enddo

deallocate(key_ids)

nfiles = get_num_blocks(diag_yaml_id, "diag_files")
allocate(file_ids(nfiles))
call get_block_ids(diag_yaml_id, "diag_files", file_ids)
print *, ""

do i = 1, nfiles
   print *, "File number:", i

   nkeys = get_nkeys(diag_yaml_id, file_ids(i))
   allocate(key_ids(nkeys))
   call get_key_ids(diag_yaml_id, file_ids(i), key_ids)

   do j = 1, nkeys
      call get_key_name(diag_yaml_id, key_ids(j), key_name)
      call get_key_value(diag_yaml_id, key_ids(j), key_value)
      print *, "  Key:", trim(key_name), " Value:", trim(key_value)
   enddo

   deallocate(key_ids)
   print *, ""
   !< The number of variables that are part of the current file
   nvariables = get_num_blocks(diag_yaml_id, "varlist", parent_block_id=file_ids(i))
   allocate(var_ids(nvariables))
   call get_block_ids(diag_yaml_id, "varlist", var_ids, parent_block_id=file_ids(i))

   do j = 1, nvariables
       print *, "  Variable number:", j

       nkeys = get_nkeys(diag_yaml_id, var_ids(j))
       allocate(key_ids(nkeys))
       call get_key_ids(diag_yaml_id, var_ids(j), key_ids)

       do k = 1, nkeys
          call get_key_name(diag_yaml_id, key_ids(k), key_name)
          call get_key_value(diag_yaml_id, key_ids(k), key_value)
          print *, "     Key:", trim(key_name), " Value:", trim(key_value)
       enddo

       deallocate(key_ids)
       print *, ""
   end do

   deallocate(var_ids)
   print *, ""
enddo
deallocate(file_ids)
call fms_end

#endif

end program parser_demo
