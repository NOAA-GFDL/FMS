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

!> @brief This program tests the diag_yaml_object_init and diag_yaml_object_end subroutines
!! in fms_diag_yaml_mod
program test_diag_yaml

#ifdef use_yaml
use FMS_mod, only: fms_init, fms_end
use fms_diag_yaml_mod
use fms_diag_yaml_object_mod
use mpp_mod
use platform_mod

implicit none

!< @brief Interface used to compare two different values
interface compare_result
subroutine compare_result_0d(key_name, res, expected_res)
  character(len=*), intent(in) :: key_name  !< Name of the key to compare
  class(*), intent(in) :: res               !< Value obtained from reading the file
  class(*), intent(in) :: expected_res      !< Value expected
end subroutine compare_result_0d

subroutine compare_result_1d(key_name, res, expected_res)
  character(len=*), intent(in) :: key_name  !< Name of the key to compare
  class(*), intent(in) :: res(:)            !< Value obtained from reading the file
  class(*), intent(in) :: expected_res(:)   !< Value expected
end subroutine compare_result_1d
end interface compare_result

type(diagYamlObject_type) :: my_yaml !< diagYamlObject obtained from diag_yaml_object_init
type(diagYamlObject_type) :: ans     !< expected diagYamlObject
logical :: checking_crashes = .false.!< Flag indicating that you are checking crashes
integer :: i !< For do loops
integer :: io_status !< The status after reading the input.nml

namelist / check_crashes_nml / checking_crashes

call fms_init()

read (input_nml_file, check_crashes_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>check_crashes: Error reading input.nml')

call diag_yaml_object_init

my_yaml = get_diag_yaml_obj()
ans = fill_ans()

if (.not. checking_crashes) then
  call compare_result("base_date", my_yaml%get_basedate(), ans%get_basedate())
  call compare_result("title", my_yaml%get_title(), ans%get_title())

  call compare_result("nfiles", size(my_yaml%diag_files, 1), size(ans%diag_files, 1))
  do i = 1, size(my_yaml%diag_files, 1)
    print *, "Working on file:", i
    call compare_diag_files(my_yaml%diag_files(i), ans%diag_files(i))
  enddo

  call compare_result("nvariables", size(my_yaml%diag_fields, 1), size(ans%diag_fields, 1))
  do i = 1, size(my_yaml%diag_fields, 1)
    print *, "Working on variable:", i
    call compare_diag_fields(my_yaml%diag_fields(i), ans%diag_fields(i))
  enddo
endif
call diag_yaml_object_end

!< Check if everything is destroyed correctly!
my_yaml = get_diag_yaml_obj()
if (allocated(my_yaml%diag_files)) call mpp_error(FATAL, "diag_file is still allocated")
if (allocated(my_yaml%diag_fields)) call mpp_error(FATAL, "diag_fields is still allocated")
call fms_end()

contains

!> @brief Compares a diagYamlFilesVar_type with the expected result and
!! crashes if they don't match
subroutine compare_diag_fields(res, expected_res)
  type(diagYamlFilesVar_type), intent(in) :: res          !< diag_field info read from yaml file
  type(diagYamlFilesVar_type), intent(in) :: expected_res !< diag_field info expected

  integer :: i !< For do loops

  call compare_result("var_fname", res%get_var_fname(), expected_res%get_var_fname())
  call compare_result("var_varname", res%get_var_varname(), expected_res%get_var_varname())
  call compare_result("var_reduction", res%get_var_reduction(), expected_res%get_var_reduction())
  call compare_result("var_module", res%get_var_module(), expected_res%get_var_module())
  call compare_result("var_skind", res%get_var_skind(), expected_res%get_var_skind())
  call compare_result("var_write", res%get_var_write(), expected_res%get_var_write())
  call compare_result("var_outname", res%get_var_outname(), expected_res%get_var_outname())
  call compare_result("var_longname", res%get_var_longname(), expected_res%get_var_longname())
!< VAR_UNITS

  call compare_result("number_variable_attributes", size(res%var_attributes,1), size(expected_res%var_attributes,1))
  do i = 1, size(res%var_attributes,1)
    call compare_attributes(res%var_attributes(i,:), expected_res%var_attributes(i,:))
  enddo
end subroutine

!> @brief Compares a diagYamlFiles_type with the expected result and
!! crashes if they don't match
subroutine compare_diag_files(res, expected_res)
  type(diagYamlFiles_type), intent(in) :: res           !< diag_file info read from yaml file
  type(diagYamlFiles_type), intent(in) :: expected_res  !< diag_file info expected

  integer :: i !< For Do loops

  call compare_result("file_fname", res%get_file_fname(), expected_res%get_file_fname())
  call compare_result("file_frequnit", res%get_file_frequnit(), expected_res%get_file_frequnit())
  call compare_result("file_freq", res%get_file_freq(), expected_res%get_file_freq())
  call compare_result("file_timeunit", res%get_file_timeunit(), expected_res%get_file_timeunit())
  call compare_result("file_unlimdim", res%get_file_unlimdim(), expected_res%get_file_unlimdim())
  call compare_result("file_write", res%get_file_write(), expected_res%get_file_write())
  call compare_result("file_realm", res%get_file_realm(), expected_res%get_file_realm())

  call compare_result("sub_region_type", res%file_sub_region%grid_type, expected_res%file_sub_region%grid_type)
  if (trim(res%file_sub_region%grid_type) .eq. "latlon") then
    call compare_result("sub_region_lat_lon", res%file_sub_region%lat_lon_sub_region, &
      expected_res%file_sub_region%lat_lon_sub_region)
  elseif (trim(res%file_sub_region%grid_type) .eq. "index") then
    call compare_result("sub_region_index", res%file_sub_region%index_sub_region, &
      expected_res%file_sub_region%index_sub_region)
    call compare_result("tile", res%file_sub_region%tile, expected_res%file_sub_region%tile)
  endif

  call compare_result("file_new_file_freq", res%get_file_new_file_freq(), expected_res%get_file_new_file_freq())
  call compare_result("file_new_file_freq_units", res%get_file_new_file_freq_units(), &
    expected_res%get_file_new_file_freq_units())
  call compare_result("file_start_time", res%get_file_start_time(), expected_res%get_file_start_time())
  call compare_result("file_duration", res%get_file_duration(), expected_res%get_file_duration())
  call compare_result("file_duration_units", res%get_file_duration_units(), expected_res%get_file_duration_units())

  call compare_result("number_global_attributes", size(res%file_global_meta,1), size(expected_res%file_global_meta,1))
  do i = 1, size(res%file_global_meta,1)
    call compare_attributes(res%file_global_meta(i,:), expected_res%file_global_meta(i,:))
  enddo

  !< Check file_varlist
  call compare_result("number_file_variables", size(res%file_varlist,1), size(expected_res%file_varlist,1))
  do i = 1, size(res%file_varlist,1)
    call compare_result("file_varlist", res%file_varlist(i), expected_res%file_varlist(i))
  enddo

end subroutine compare_diag_files

!< @brief Compares a key/value pair of attributes (global or variable)
subroutine compare_attributes(res, expected_res)
  character (len=255), dimension(:), intent(in) :: res           !< res(1) key (attribute name) read from file
                                                                 !< res(2) value (attribute value) read from file
  character (len=255), dimension(:), intent(in) :: expected_res  !< res(1) key (attribute name) expected
                                                                 !< res(2) value (attribute value) expected

  call compare_result("attributes key", res(1), expected_res(1))
  call compare_result("attributes value", res(2), expected_res(2))

end subroutine compare_attributes

!< @brief Fill in a diag_yaml_object with the expected result
function fill_ans() &
result(ans)
  type(diagYamlObject_type) :: ans

  ans%diag_basedate = (/2, 1, 1, 0, 0 , 0 /)
  ans%diag_title = "test_diag_manager"

  !< CHECK DIAG_FILES
  allocate(ans%diag_files(3))
  ans%diag_files(1)%file_fname = "wild_card_name%4yr%2mo%2dy%2hr"
  ans%diag_files(2)%file_fname = "normal"
  ans%diag_files(3)%file_fname = "normal2"

  ans%diag_files(1)%file_freq = 6
  ans%diag_files(2)%file_freq = 24
  ans%diag_files(3)%file_freq = 24

  ans%diag_files(1)%file_frequnit = "hours"
  ans%diag_files(2)%file_frequnit = "days"
  ans%diag_files(3)%file_frequnit = "days"

  ans%diag_files(1)%file_timeunit = "hours"
  ans%diag_files(2)%file_timeunit = "hours"
  ans%diag_files(3)%file_timeunit = "hours"

  ans%diag_files(1)%file_unlimdim = "time"
  ans%diag_files(2)%file_unlimdim = "records"
  ans%diag_files(3)%file_unlimdim = "records"

  !< Check if the realm is set correctly
  ans%diag_files(1)%file_realm = "ATM"

  !< Check if file_write is set correctly
  ans%diag_files(1)%file_write = .false.
  !< The second file does not have a "write_file" entry, and the default is true!
  ans%diag_files(2)%file_write = .true.
  !< The third file has the "write_file" entry set as true
  ans%diag_files(3)%file_write = .true.

  !< Check if the latlon region is set correctly
  ans%diag_files(2)%file_sub_region%grid_type = "latlon"
  allocate(ans%diag_files(2)%file_sub_region%lat_lon_sub_region(8))
  ans%diag_files(2)%file_sub_region%lat_lon_sub_region = (/64., -999., -999., -999., -999., 20., -999., -999. /)

  !< Check if the index region is set correctly
  ans%diag_files(3)%file_sub_region%grid_type = "index"
  allocate(ans%diag_files(3)%file_sub_region%index_sub_region(8))
  ans%diag_files(3)%file_sub_region%index_sub_region = (/10, -999, 10, 20, -999, -999, -999, -999 /)
  ans%diag_files(3)%file_sub_region%tile = 1

  !< Check if the file_new_* stuff is set correctly
  ans%diag_files(1)%file_new_file_freq = 6
  ans%diag_files(1)%file_new_file_freq_units = "hours"
  ans%diag_files(1)%file_start_time = "2 1 1 0 0 0"

  !< Check the file_duration stuff is set correctly
  ans%diag_files(1)%file_duration = 12
  ans%diag_files(1)%file_duration_units = "hours"

  !< The second file is not using the new_* stuff, so it should be zero
  ans%diag_files(2)%file_new_file_freq = 0
  ans%diag_files(2)%file_duration = 0

  !< Check if the global attribute block is read correctly
  allocate(ans%diag_files(1)%file_global_meta(1, 2))
  ans%diag_files(1)%file_global_meta(1, 1) = "is_a_file"
  ans%diag_files(1)%file_global_meta(1, 2) = "true"

  !< Check if the file varlist is set correctly!
  allocate(ans%diag_files(1)%file_varlist(1))
  ans%diag_files(1)%file_varlist(1) = "sst"
  allocate(ans%diag_files(2)%file_varlist(1))
  ans%diag_files(2)%file_varlist(1) = "sst"
  allocate(ans%diag_files(3)%file_varlist(1))
  ans%diag_files(3)%file_varlist(1) = "sst"

  !< CHECK DIAG_FIELDS
  allocate(ans%diag_fields(3))
  do i = 1, size(ans%diag_fields)
    call fill_field_defaults(ans%diag_fields(i))
  enddo

  ans%diag_fields(1)%var_fname="wild_card_name%4yr%2mo%2dy%2hr"
  ans%diag_fields(2)%var_fname="normal"
  ans%diag_fields(3)%var_fname="normal2"

  !< Test if the attribute block was read correctly:
  allocate(ans%diag_fields(2)%var_attributes(1, 2))
  ans%diag_fields(2)%var_attributes(1, 1) = "do_sst"
  ans%diag_fields(2)%var_attributes(1, 2) = ".true."

  !< Test if the var_write was set correctly:
  ans%diag_fields(1)%var_write = .false.
  ans%diag_fields(2)%var_write = .true.
  ans%diag_fields(3)%var_write = .true. !< Here the write_var key does not exist and the default should be true

  !< Test if the optional var_longname was read correctly:
  ans%diag_fields(3)%var_longname = "S S T"

end function fill_ans

!< @brief Fills in diagYamlFilesVar_type with the expected results
subroutine fill_field_defaults(ans)
  type(diagYamlFilesVar_type) :: ans !< diag_field type to fill

  ans%var_varname = "sst"
  ans%var_reduction = "average"
  ans%var_module = "test_diag_manager_mod"
  ans%var_skind = "float"
  ans%var_outname = "sst"
  ans%var_write = .true.

end subroutine fill_field_defaults
#endif
end program test_diag_yaml

#ifdef use_yaml
!< @brief Compare a key value with the expected result
subroutine compare_result_0d(key_name, res, expected_res)
  use platform_mod
  use mpp_mod
  character(len=*), intent(in) :: key_name  !< Name of the key to compare
  class(*), intent(in) :: res               !< Value obtained from reading the file
  class(*), intent(in) :: expected_res      !< Value expected

  print *, "Comparing ", trim(key_name)
  select type(res)
    type is(character(len=*))
      select type(expected_res)
        type is(character(len=*))
          if(trim(res) .ne. trim(expected_res)) &
            call mpp_error(FATAL, "Error!: "//trim(key_name)//" is not the expected result. "//trim(res)//" ne "//&
              trim(expected_res)//".")
      end select
    type is (integer(kind=i4_kind))
       select type(expected_res)
         type is(integer(kind=i4_kind))
           if (res .ne. expected_res) then
             print *, res, " ne ", expected_res
            call mpp_error(FATAL, "Error!: "//trim(key_name)//" is not the expected result.")
          endif
      end select
    type is (logical)
       select type(expected_res)
         type is(logical)
           if ((res .and. .not. expected_res) .or. (.not. res .and. expected_res)) then
            print*, res, " ne ", expected_res
            call mpp_error(FATAL, "Error!:"//trim(key_name)//" is not the expected result")
           endif
      end select
  end select

end subroutine compare_result_0d

!< @brief Compare a 1d key value with the expected result
subroutine compare_result_1d(key_name, res, expected_res)
  use platform_mod
  use mpp_mod
  character(len=*), intent(in) :: key_name  !< Name of the key to compare
  class(*), intent(in) :: res(:)            !< Value obtained from reading the file
  class(*), intent(in) :: expected_res(:)   !< Value expected

  integer :: i

  print *, "Comparing ", trim(key_name)

  select type(res)
    type is (integer(kind=i4_kind))
      select type(expected_res)
        type is (integer(kind=i4_kind))
          do i = 1, size(res,1)
            if( res(i) .ne. expected_res(i)) then
              print *, res, " ne ", expected_res
              call mpp_error(FATAL, "Error!: "//trim(key_name)//" is not the expected result. ")
            endif
          enddo
      end select
    type is (real(kind=r4_kind))
      select type(expected_res)
        type is (real(kind=r4_kind))
          do i = 1, size(res,1)
            if( res(i) .ne. expected_res(i)) then
              print *, res, " ne ", expected_res
              call mpp_error(FATAL, "Error!: "//trim(key_name)//" is not the expected result. ")
            endif
          enddo
      end select
    type is (real(kind=r8_kind))
      select type(expected_res)
        type is (real(kind=r8_kind))
          do i = 1, size(res,1)
            if( res(i) .ne. expected_res(i)) then
              print *, res, " ne ", expected_res
              call mpp_error(FATAL, "Error!: "//trim(key_name)//" is not the expected result. ")
            endif
          enddo
    end select
  end select
end subroutine compare_result_1d
#endif