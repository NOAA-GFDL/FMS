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

!> @brief This program tests the diag_model_subset feature of diag_mananger_init
!! It requires two PEs to run and it runs with diag_table_yaml_27
program test_diag_ocean

#ifdef use_yaml
use FMS_mod, only: fms_init, fms_end, string
use fms_diag_yaml_mod
use diag_manager_mod, only: diag_manager_init
use diag_data_mod, only: DIAG_NULL, DIAG_OCEAN, DIAG_OTHER
use mpp_mod
use platform_mod

implicit none

type(diagYamlObject_type) :: my_yaml !< diagYamlObject obtained from diag_yaml_object_init
type(diagYamlFiles_type), allocatable, dimension (:) :: diag_files !< Files from the diag_yaml
type(diagYamlFilesVar_type), allocatable, dimension(:) :: diag_fields !< Fields from the diag_yaml
character(len=10), allocatable :: file_names(:) !< The expected names of the files
character(len=10), allocatable :: var_names(:) !< The expected names of the variables
integer :: diag_subset !< Diag_subset to be sent to diag_manager_init
integer :: nfiles !< Expected number of files
integer :: nvariables !< Expected number of variables
integer :: i !< For do loops

call fms_init()

if (mpp_npes() .ne. 2) call mpp_error(FATAL, "test_diag_ocean requires two PEs!")

!> PE 0 is not going to include the file with is_ocean = .true.
if (mpp_pe() .eq. 0) then
  diag_subset = DIAG_OTHER
  nfiles = 2
  allocate(file_names(nfiles))
  file_names = (/"file1", "file3"/)
  nvariables = 3
  allocate(var_names(nvariables))
  var_names = (/"sst1", "sst3", "sst4"/)
endif

!> PE 1 is only going to include the file with is_ocean = .true.
if (mpp_pe() .eq. 1) then
  diag_subset = DIAG_OCEAN
  nfiles = 1
  allocate(file_names(nfiles))
  file_names = (/"file2"/)
  nvariables = 1
  allocate(var_names(nvariables))
  var_names = (/"sst2"/)
endif

call diag_manager_init(diag_model_subset=diag_subset)

my_yaml = get_diag_yaml_obj()
diag_files = my_yaml%get_diag_files()
if (size(diag_files) .ne. nfiles) call mpp_error(FATAL, "The number of files should be "//string(nfiles))

do i = 1, nfiles
  if(trim(file_names(i)) .ne. diag_files(i)%get_file_fname()) &
    call mpp_error(FATAL, "The file_name should of the "//string(i)//" file should be "//&
      &trim(file_names(i))//" not "//diag_files(i)%get_file_fname())
end do

diag_fields = my_yaml%get_diag_fields()
if (size(diag_fields) .ne. nvariables) call mpp_error(FATAL, "The number of variables should be "//string(nvariables))

do i = 1, nvariables
  if(trim(var_names(i)) .ne. diag_fields(i)%get_var_varname()) &
    call mpp_error(FATAL, "The var_name should of the "//string(i)//" field should be "//&
      &trim(var_names(i))//" not "//diag_fields(i)%get_var_varname())
end do

deallocate(diag_files)
deallocate(diag_fields)
deallocate(file_names)
deallocate(var_names)

call diag_yaml_object_end
call fms_end()

#endif
end program test_diag_ocean