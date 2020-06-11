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

!> @file
!! @brief Unit test for the STDIN function
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_stdin
  use mpp_mod, only : mpp_init, mpp_exit
  use mpp_mod, only : mpp_error, FATAL
  use mpp_mod, only : stdin
  use iso_fortran_env, only : INPUT_UNIT

  integer :: in_unit !< Stores the returned standard input unit number
  character(len=128) :: err_msg !< Stores an arbitrary error message
  logical :: open5 !< Indicates whether unit 5 is open
  logical :: open100 !< Indicates whether unit 100 is open

  in_unit = stdin()

  write(*,*) "stdin() should get and return the value of in_unit."
  write(*,*) "This value should match INPUT_UNIT from iso_fortran_env."

  if (INPUT_UNIT.eq.in_unit) then
    write(*,*) "PASS: stdin() returned the correct value of in_unit."
  else
    write(*,*) "Integer returned by stdin():"
    write(*,*) in_unit
    write(*,*) "Integer expected:"
    write(*,*) INPUT_UNIT
    call mpp_error(FATAL, "ERROR: stdin() did not return the expected value")
  end if

end program test_stdin
