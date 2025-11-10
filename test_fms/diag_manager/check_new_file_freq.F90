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

!> @brief  General program to check the output corner cases with the new_file_freq keys
program test_check_new_file_freq
  use fms_mod, only: fms_init, fms_end
  use fms2_io_mod,      only: FmsNetcdfFile_t, open_file, close_file, get_dimension_size
  use mpp_mod,          only: FATAL, mpp_error
  implicit none

  type(FmsNetcdfFile_t)              :: fileobj     !< FMS2 fileobj
  integer                            :: var_size    !< Size of the time dimension

  call fms_init()
  ! Check the output generated from this diag_table yaml
  ! title: test_none
  ! base_date: 2 1 1 0 0 0
  ! diag_files:
  ! - file_name: test_none%4yr%2mo%2dy%2hr
  !   freq: 6 hours
  !   time_units: hours
  !   unlimdim: time
  !   new_file_freq: 12 hours
  !   file_duration: 18 hours
  !   varlist:
  !   - module: ocn_mod
  !     var_name: var0
  !     output_name: var0_none
  !     reduction: none
  !     kind: r4

  if (.not. open_file(fileobj, "test_none_0002_01_01_06.nc", "read")) &
    call mpp_error(FATAL, "unable to open test_none_0002_01_01_06.nc")

  call get_dimension_size(fileobj, "time", var_size)
  if (var_size .ne. 2) call mpp_error(FATAL, "The dimension of time in the file:test_none_0002_01_01_06.nc"//&
                                                 " is not the correct size!")

  call close_file(fileobj)

  if (.not. open_file(fileobj, "test_none_0002_01_01_15.nc", "read")) &
    call mpp_error(FATAL, "unable to open test_none_0002_01_01_15.nc")

  call get_dimension_size(fileobj, "time", var_size)
  if (var_size .ne. 1) call mpp_error(FATAL, "The dimension of time in the file:test_none_0002_01_01_15.nc"//&
                                                 " is not the correct size!")

  call close_file(fileobj)

  call fms_end()
end program test_check_new_file_freq