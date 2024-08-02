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

!> @brief  This programs tests the subroutines get_num_unique_blocks, get_unique_block_ids, and
!! get_block_name
program generic_blocks
#ifdef use_yaml
  use fms_mod, only: fms_init, fms_end
  use mpp_mod, only: mpp_error, FATAL
  use yaml_parser_mod

  implicit none

  integer              :: yaml_id            !< Id of the yaml file
  integer, allocatable :: field_table_ids(:) !< The Ids of the field table entries
  integer, allocatable :: modlist_ids(:)     !< The ids of the mods entries
  integer, allocatable :: varlist_ids(:)     !< The ids of the variable entries
  integer, allocatable :: block_ids(:)       !< The ids of the block entries
  integer, allocatable :: misc_block_ids(:)  !< The ids of the misc block entries
  integer, allocatable :: key_ids(:)         !< The ids of the keys
  character(len=50)    :: variable_name      !< The variable name
  character(len=50)    :: model_type_name    !< The model type
  character(len=50)    :: block_name         !< The name of the block
  character(len=50)    :: key_name           !< The name of the key
  character(len=50)    :: key_value          !< The value of the key
  character(len=50)    :: varnames(2)        !< The expected names of the variables
  character(len=50)    :: blocknames1(1)     !< The expected names of the blocks for the first variable
  character(len=50)    :: blocknames2(2)     !< The expected names of the blocks for the second variable
  character(len=50)    :: keys(5)            !< The expected names of the keys
  character(len=50)    :: values(5)          !< The expected names values of they keys
  integer              :: key_count          !< To keep track of the expected answers

  logical :: correct_answer   !< True if the answer is correct
  integer :: i, j, k, l, m, n !< For do loops

  call fms_init()
  varnames(1) = "sphum"
  varnames(2) = "soa"

  blocknames1(1) = "profile_type"
  blocknames2(1) = "chem_param"
  blocknames2(2) = "profile_type"

  key_count = 0
  keys(1) = "value";         values(1) = "fixed"
  keys(2) = "surface_value"; values(2) = "3.0e-06"
  keys(3) = "value";         values(3) = "aerosol"
  keys(4) = "value";         values(4) = "fixed"
  keys(5) = "surface_value"; values(5) = "1.0e-32"

  yaml_id = open_and_parse_file("sample.yaml")
  allocate(field_table_ids(get_num_blocks(yaml_id, "field_table")))
  call get_block_ids(yaml_id, "field_table", field_table_ids)
  do i = 1, size(field_table_ids)
    allocate(modlist_ids(get_num_blocks(yaml_id, "modlist", parent_block_id=field_table_ids(i))))
    call get_block_ids(yaml_id, "modlist", modlist_ids, field_table_ids(i))

    do j = 1, size(modlist_ids)
      call get_value_from_key(yaml_id, modlist_ids(j), "model_type", model_type_name)
      print *, "Modlist::", trim(model_type_name)
      if (trim(model_type_name) .ne. "atmos_mod") &
        call mpp_error(FATAL, "Modlist is not the expected result")

      allocate(varlist_ids(get_num_blocks(yaml_id, "varlist", parent_block_id=modlist_ids(j))))
      call get_block_ids(yaml_id, "varlist", varlist_ids, modlist_ids(j))

      do k = 1, size(varlist_ids)
        call get_value_from_key(yaml_id, varlist_ids(k), "variable", variable_name)
        print *, "Variable::", trim(variable_name)
        if (trim(variable_name) .ne. varnames(k)) &
          call mpp_error(FATAL, "Variable is not the expected result")

        allocate(block_ids(get_num_unique_blocks(yaml_id, parent_block_id=varlist_ids(k))))
        call get_unique_block_ids(yaml_id, block_ids, parent_block_id=varlist_ids(k))
        do l = 1, size(block_ids)
          call get_block_name(yaml_id, block_ids(l), block_name)
          print *, "Block_name::", trim(block_name)

          if (k == 1) then
            correct_answer = trim(blocknames1(l)) .eq. trim(block_name)
          else
            correct_answer = trim(blocknames2(l)) .eq. trim(block_name)
          endif

          if (.not. correct_answer) call mpp_error(FATAL, "blockname is not the expected result")
          allocate(misc_block_ids(get_num_blocks(yaml_id, block_name, parent_block_id=varlist_ids(k))))
          call get_block_ids(yaml_id, block_name, misc_block_ids, parent_block_id=varlist_ids(k))
          do m = 1, size(misc_block_ids)
            allocate(key_ids(get_nkeys(yaml_id, misc_block_ids(m))))
            call get_key_ids(yaml_id, misc_block_ids(m), key_ids)
            do n = 1, size(key_ids)
              key_count = key_count + 1
              call get_key_name(yaml_id, key_ids(n), key_name)
              call get_key_value(yaml_id, key_ids(n), key_value)
              print *, "KEY:", trim(key_name), " VALUE:", trim(key_value)

              if (trim(key_name) .ne. trim(keys(key_count))) &
                call mpp_error(FATAL, "The key is not correct")

              if (trim(key_value) .ne. trim(values(key_count))) &
                call mpp_error(FATAL, "The value is not correct")
            enddo
            deallocate(key_ids)
          enddo
          deallocate(misc_block_ids)
        enddo
        deallocate(block_ids)
        print *, "---------"
      enddo
      deallocate(varlist_ids)
    enddo
    deallocate(modlist_ids)
  enddo
  call fms_end()
#endif
end program generic_blocks
