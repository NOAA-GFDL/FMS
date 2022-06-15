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

use FMS_mod, only: fms_init, fms_end
use fms_diag_yaml_mod
use diag_data_mod, only: DIAG_NULL, DIAG_ALL, get_base_year, get_base_month, get_base_day, get_base_hour, &
                       & get_base_minute, get_base_second, diag_data_init
use  time_manager_mod, only: set_calendar_type, JULIAN
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

logical :: checking_crashes = .false.!< Flag indicating that you are checking crashes
integer :: i !< For do loops
integer :: io_status !< The status after reading the input.nml
integer, allocatable :: indices(:) !< Array of indices

#ifdef use_yaml
type(diagYamlFiles_type), allocatable, dimension (:) :: diag_files !< Files from the diag_yaml
type(diagYamlFilesVar_type), allocatable, dimension(:) :: diag_fields !< Fields from the diag_yaml
type(diagYamlObject_type) :: my_yaml !< diagYamlObject obtained from diag_yaml_object_init
type(diagYamlObject_type) :: ans     !< expected diagYamlObject
#endif

namelist / check_crashes_nml / checking_crashes

call fms_init()

read (input_nml_file, check_crashes_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>check_crashes: Error reading input.nml')

#ifndef use_yaml
if (checking_crashes) call mpp_error(FATAL, "It is crashing!")
call fms_end()
#else

call set_calendar_type(JULIAN)
call diag_data_init()
call diag_yaml_object_init(DIAG_ALL)

my_yaml = get_diag_yaml_obj()

if (.not. checking_crashes) then
  call compare_result("base_date", my_yaml%get_basedate(), (/2, 1, 1, 0, 0 , 0 /))
  call check_base_time()

  call compare_result("title", my_yaml%get_title(), "test_diag_manager")

  diag_files = my_yaml%get_diag_files()
  call compare_result("nfiles", size(diag_files), 3) !< the fourth file has file_write = false so it doesn't count
  call compare_diag_files(diag_files)

  diag_fields = my_yaml%get_diag_fields()
  call compare_result("nfields", size(diag_fields), 3) !< the fourth variable has var_write = false so it doesn't count
  call compare_diag_fields(diag_fields)

  !< Check that get_num_unique_fields is getting the correct number of unique fields
  call compare_result("number of unique fields", get_num_unique_fields(), 2)

  deallocate(diag_files)
  deallocate(diag_fields)

  indices = find_diag_field("sst")
  print *, "sst was found in ", indices
  if (size(indices) .ne. 2) &
    call mpp_error(FATAL, 'sst was supposed to be found twice!')
  if (indices(1) .ne. 2 .and. indices(2) .ne. 1) &
    call mpp_error(FATAL, 'sst was supposed to be found in indices 1 and 2')

  diag_fields = get_diag_fields_entries(indices)
  call compare_result("sst - nfields", size(diag_fields), 2)
  call compare_result("sst - fieldname", diag_fields(1)%get_var_varname(), "sst")
  call compare_result("sst - fieldname", diag_fields(2)%get_var_varname(), "sst")
  deallocate(diag_fields)

  diag_files = get_diag_files_entries(indices)
  call compare_result("sst - nfiles", size(diag_files), 2)
  call compare_result("sst - filename", diag_files(1)%get_file_fname(), "normal")
  call compare_result("sst - filename", diag_files(2)%get_file_fname(), "wild_card_name%4yr%2mo%2dy%2hr")
  deallocate(diag_files)
  deallocate(indices)

  indices = find_diag_field("sstt")
  print *, "sstt was found in ", indices
  if (size(indices) .ne. 1) &
    call mpp_error(FATAL, 'sstt was supposed to be found twice!')
  if (indices(1) .ne. 3) &
    call mpp_error(FATAL, 'sstt was supposed to be found in indices 1 and 2')
  deallocate(indices)

  indices = find_diag_field("sstt2") !< This is in diag_table but it has write_var = false
  print *, "sstt2 was found in ", indices
  if (indices(1) .ne. -999) &
    call mpp_error(FATAL, "sstt2 is not in the diag_table!")

  indices = find_diag_field("tamales")
  print *, "tamales was found in ", indices
  if (indices(1) .ne. -999) &
    call mpp_error(FATAL, "tamales is not in the diag_table!")

endif

call diag_yaml_object_end

call fms_end()

contains

!> @brief Compares a diagYamlFilesVar_type with the expected result and
!! crashes if they don't match
subroutine compare_diag_fields(res)
  type(diagYamlFilesVar_type), intent(in) :: res(:)          !< diag_field info read from yaml file
  character (len=255), dimension(:, :), allocatable :: var_attributes !< Variable attributes

  call compare_result("var_fname 1", res(1)%get_var_fname(), "wild_card_name%4yr%2mo%2dy%2hr")
  call compare_result("var_fname 2", res(2)%get_var_fname(), "normal")
  call compare_result("var_fname 3", res(3)%get_var_fname(), "normal2")

  call compare_result("var_varname 1", res(1)%get_var_varname(), "sst")
  call compare_result("var_varname 2", res(2)%get_var_varname(), "sst")
  call compare_result("var_varname 3", res(3)%get_var_varname(), "sstt")

  call compare_result("var_reduction 1", res(1)%get_var_reduction(), "average")
  call compare_result("var_reduction 2", res(2)%get_var_reduction(), "average")
  call compare_result("var_reduction 3", res(3)%get_var_reduction(), "average")

  call compare_result("var_module 1", res(1)%get_var_module(), "test_diag_manager_mod")
  call compare_result("var_module 2", res(2)%get_var_module(), "test_diag_manager_mod")
  call compare_result("var_module 3", res(3)%get_var_module(), "test_diag_manager_mod")

  call compare_result("var_skind 1", res(1)%get_var_skind(), "r4")
  call compare_result("var_skind 2", res(2)%get_var_skind(), "r4")
  call compare_result("var_skind 3", res(3)%get_var_skind(), "r4")

  call compare_result("var_outname 1", res(1)%get_var_outname(), "sst")
  call compare_result("var_outname 2", res(2)%get_var_outname(), "sst")
  call compare_result("var_outname 3", res(3)%get_var_outname(), "sstt")

  call compare_result("var_longname 1", res(1)%get_var_longname(), "")
  call compare_result("var_longname 2", res(2)%get_var_longname(), "")
  call compare_result("var_longname 3", res(3)%get_var_longname(), "S S T")

  if (res(1)%is_var_attributes()) call mpp_error(FATAL, "The variable attributes for the first file was set?")

  var_attributes = res(2)%get_var_attributes()
  if (.not. allocated(var_attributes)) call mpp_error(FATAL, "The variable attributes for the second file was not set")
  call compare_result("var attributes key", var_attributes(1,1), "do_sst")
  call compare_result("var attributes value", var_attributes(1,2), ".true.")
  deallocate(var_attributes)

  if (res(3)%is_var_attributes()) call mpp_error(FATAL, "The variable attributes for the third file was set?")

end subroutine

!> @brief Compares a diagYamlFiles_type with the expected result and
!! crashes if they don't match
subroutine compare_diag_files(res)
  type(diagYamlFiles_type), intent(in) :: res(:)           !< diag_file info read from yaml file

  character (len=255), dimension(:), allocatable :: varlist !< List of variables
  character (len=255), dimension(:, :), allocatable :: global_meta !< List of global meta

  call compare_result("file_fname 1", res(1)%get_file_fname(), "wild_card_name%4yr%2mo%2dy%2hr")
  call compare_result("file_fname 2", res(2)%get_file_fname(), "normal")
  call compare_result("file_fname 3", res(3)%get_file_fname(), "normal2")

  call compare_result("file_freq 1", res(1)%get_file_freq(), 6)
  call compare_result("file_freq 2", res(2)%get_file_freq(), 24)
  call compare_result("file_freq 3", res(3)%get_file_freq(), -1)

  call compare_result("file_frequnit 1", res(1)%get_file_frequnit(), "hours")
  call compare_result("file_frequnit 2", res(2)%get_file_frequnit(), "days")
  call compare_result("file_frequnit 3", res(3)%get_file_frequnit(), "days")

  call compare_result("file_timeunit 1", res(1)%get_file_timeunit(), "hours")
  call compare_result("file_timeunit 2", res(2)%get_file_timeunit(), "hours")
  call compare_result("file_timeunit 3", res(3)%get_file_timeunit(), "hours")

  call compare_result("file_unlimdim 1", res(1)%get_file_unlimdim(), "time")
  call compare_result("file_unlimdim 2", res(2)%get_file_unlimdim(), "records")
  call compare_result("file_unlimdim 3", res(3)%get_file_unlimdim(), "records")

  call compare_result("file_new_file_freq 1", res(1)%get_file_new_file_freq(), 6)
  call compare_result("file_new_file_freq 2", res(2)%get_file_new_file_freq(), DIAG_NULL)
  call compare_result("file_new_file_freq 3", res(3)%get_file_new_file_freq(), DIAG_NULL)

  call compare_result("file_new_file_freq_units 1", res(1)%get_file_new_file_freq_units(), "hours")
  call compare_result("file_new_file_freq_units 2", res(2)%get_file_new_file_freq_units(), "")
  call compare_result("file_new_file_freq_units 3", res(3)%get_file_new_file_freq_units(), "")

  call compare_result("file_duration 1", res(1)%get_file_duration(), 12)
  call compare_result("file_duration 2", res(2)%get_file_duration(), DIAG_NULL)
  call compare_result("file_duration 3", res(3)%get_file_duration(), DIAG_NULL)

  call compare_result("file_duration_units 1", res(1)%get_file_duration_units(), "hours")
  call compare_result("file_duration_units 2", res(2)%get_file_duration_units(), "")
  call compare_result("file_duration_units 3", res(3)%get_file_duration_units(), "")

  call compare_result("file_start_time 1", res(1)%get_file_start_time(), "2 1 1 0 0 0")
  call compare_result("file_start_time 2", res(2)%get_file_start_time(), "")
  call compare_result("file_start_time 3", res(3)%get_file_start_time(), "")

  varlist = res(1)%get_file_varlist()
  if (.not. allocated(varlist)) call mpp_error(FATAL, "The varlist for the first file was not set")
  call compare_result("number_variables 1", size(varlist), 1)
  call compare_result("varlist 1", varlist(1), "sst")
  deallocate(varlist)

  varlist = res(2)%get_file_varlist()
  if (.not. allocated(varlist)) call mpp_error(FATAL, "The varlist for the first file was not set")
  call compare_result("number_variables 2", size(varlist), 1)
  call compare_result("varlist 2", varlist(1), "sst")
  deallocate(varlist)

  varlist = res(3)%get_file_varlist()
  if (.not. allocated(varlist)) call mpp_error(FATAL, "The varlist for the first file was not set")
  call compare_result("number_variables 3", size(varlist), 1)
  call compare_result("varlist 3", varlist(1), "sstt")
  deallocate(varlist)

  global_meta= res(1)%get_file_global_meta()
  if (.not. allocated(global_meta)) call mpp_error(FATAL, "The global meta for the first file was not set")
  call compare_result("attributes key", global_meta(1,1), "is_a_file")
  call compare_result("attributes value", global_meta(1,2), "true")
  deallocate(global_meta)

  if (res(2)%is_global_meta()) call mpp_error(FATAL, "The global meta for the second file was set?")
  if (res(3)%is_global_meta()) call mpp_error(FATAL, "The global meta for the third file was set?")

end subroutine compare_diag_files

!> @brief Check if the base_time saved in diag_data is correct
subroutine check_base_time()
  integer :: base_time_mod_var(6) !< The base_time obtained from diag_data

  base_time_mod_var(1) = get_base_year()
  base_time_mod_var(2) = get_base_month()
  base_time_mod_var(3) = get_base_day()
  base_time_mod_var(4) = get_base_hour()
  base_time_mod_var(5) = get_base_minute()
  base_time_mod_var(6) = get_base_second()

  call compare_result("base_time", base_time_mod_var, (/2, 1, 1, 0, 0 ,0 /))
end subroutine check_base_time

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
