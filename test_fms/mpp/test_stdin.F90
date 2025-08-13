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
!! @brief Unit test for the STDIN function
!! @author Colin Gladue
program test_stdin
  use mpp_mod, only : stdin
  use iso_fortran_env, only : INPUT_UNIT, OUTPUT_UNIT

  integer :: in_unit !< Stores the returned standard input unit number

  in_unit = stdin()

  write(OUTPUT_UNIT,*) "stdin() should get and return the value of in_unit."
  write(OUTPUT_UNIT,*) "This value should match INPUT_UNIT from iso_fortran_env."

  if (INPUT_UNIT.eq.in_unit) then
    write(OUTPUT_UNIT,*) "PASS: stdin() returned the correct value of in_unit."
  else
    write(OUTPUT_UNIT,*) "Integer returned by stdin():"
    write(OUTPUT_UNIT,*) in_unit
    write(OUTPUT_UNIT,*) "Integer expected:"
    write(OUTPUT_UNIT,*) INPUT_UNIT
    write(OUTPUT_UNIT,*) "ERROR: stdin() did not return the expected value"
    stop 1
  end if

end program test_stdin
