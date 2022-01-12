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

type(diagYamlFiles_type), allocatable, dimension (:) :: diag_files !< Files from the diag_yaml
type(diagYamlFilesVar_type), allocatable, dimension(:) :: diag_fields !< Fields from the diag_yaml

namelist / check_crashes_nml / checking_crashes

call fms_init()

read (input_nml_file, check_crashes_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>check_crashes: Error reading input.nml')

call diag_yaml_object_init

my_yaml = get_diag_yaml_obj()

if (.not. checking_crashes) then
  call compare_result("base_date", my_yaml%get_basedate(), (/2, 1, 1, 0, 0 , 0 /))
  call compare_result("title", my_yaml%get_title(), "test_diag_manager")

  diag_files = my_yaml%get_diag_files()
  call compare_result("nfiles", size(diag_files), 3)
  call compare_diag_files(diag_files)

  diag_fields = my_yaml%get_diag_fields()
  call compare_result("nfields", size(diag_fields), 3)
  call compare_diag_fields(diag_fields)

endif
deallocate(diag_files)
deallocate(diag_fields)

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

  call compare_result("var_skind 1", res(1)%get_var_skind(), "float")
  call compare_result("var_skind 2", res(2)%get_var_skind(), "float")
  call compare_result("var_skind 3", res(3)%get_var_skind(), "float")

  call compare_result("var_write 1", res(1)%get_var_write(), .false.)
  call compare_result("var_write 2", res(2)%get_var_write(), .true.)
  call compare_result("var_write 3", res(3)%get_var_write(), .true.)

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
  call compare_result("file_freq 3", res(3)%get_file_freq(), 24)

  call compare_result("file_frequnit 1", res(1)%get_file_frequnit(), "hours")
  call compare_result("file_frequnit 2", res(2)%get_file_frequnit(), "days")
  call compare_result("file_frequnit 3", res(3)%get_file_frequnit(), "days")

  call compare_result("file_timeunit 1", res(1)%get_file_timeunit(), "hours")
  call compare_result("file_timeunit 2", res(2)%get_file_timeunit(), "hours")
  call compare_result("file_timeunit 3", res(3)%get_file_timeunit(), "hours")

  call compare_result("file_unlimdim 1", res(1)%get_file_unlimdim(), "time")
  call compare_result("file_unlimdim 2", res(2)%get_file_unlimdim(), "records")
  call compare_result("file_unlimdim 3", res(3)%get_file_unlimdim(), "records")

  call compare_result("file_realm 1", res(1)%get_file_realm(), "ATM")
  call compare_result("file_realm 2", res(2)%get_file_realm(), "")
  call compare_result("file_realm 3", res(3)%get_file_realm(), "")

  call compare_result("file_write 1", res(1)%get_file_write(), .false.)
  call compare_result("file_write 2", res(2)%get_file_write(), .true.)
  call compare_result("file_write 3", res(3)%get_file_write(), .true.)

  call compare_result("file_new_file_freq 1", res(1)%get_file_new_file_freq(), 6)
  call compare_result("file_new_file_freq 2", res(2)%get_file_new_file_freq(), 0)
  call compare_result("file_new_file_freq 3", res(3)%get_file_new_file_freq(), 0)

  call compare_result("file_new_file_freq_units 1", res(1)%get_file_new_file_freq_units(), "hours")
  call compare_result("file_new_file_freq_units 2", res(2)%get_file_new_file_freq_units(), "")
  call compare_result("file_new_file_freq_units 3", res(3)%get_file_new_file_freq_units(), "")

  call compare_result("file_duration 1", res(1)%get_file_duration(), 12)
  call compare_result("file_duration 2", res(2)%get_file_duration(), 0)
  call compare_result("file_duration 3", res(3)%get_file_duration(), 0)

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
