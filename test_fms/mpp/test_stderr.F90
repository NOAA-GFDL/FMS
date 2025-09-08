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

!> @file
!! @brief Unit test for the STDERR function
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_stderr
  use mpp_mod, only : stderr
  use iso_fortran_env, only : ERROR_UNIT, OUTPUT_UNIT

  integer :: err_unit !< Stores the returned standard error unit number

  err_unit = stderr()

  write(OUTPUT_UNIT,*) "stderr() should get and return the value of err_unit."
  write(OUTPUT_UNIT,*) "This value should match ERROR_UNIT from iso_fortran_env."

  if (ERROR_UNIT.eq.err_unit) then
    write(OUTPUT_UNIT,*) "PASS: stderr() returned the correct value of err_unit."
  else
    write(OUTPUT_UNIT,*) "Integer returned by stderr():"
    write(OUTPUT_UNIT,*) err_unit
    write(OUTPUT_UNIT,*) "Integer expected:"
    write(OUTPUT_UNIT,*) ERROR_UNIT
    write(OUTPUT_UNIT,*) "ERROR: stderr() did not return the expected value"
    stop 1
  end if

end program test_stderr
