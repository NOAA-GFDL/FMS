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
!! @brief Unit test for the STDERR function
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_stderr
  use mpp_mod, only : mpp_error, FATAL
  use mpp_mod, only : stderr
  use iso_fortran_env, only : ERROR_UNIT

  integer :: err_unit !< Stores the returned standard error unit number

  err_unit = stderr()

  write(*,*) "stderr() should get and return the value of err_unit."
  write(*,*) "This value should match ERROR_UNIT from iso_fortran_env."

  if (ERROR_UNIT.eq.in_unit) then
    write(*,*) "PASS: stderr() returned the correct value of err_unit."
  else
    write(*,*) "Integer returned by stderr():"
    write(*,*) err_unit
    write(*,*) "Integer expected:"
    write(*,*) ERROR_UNIT
    call mpp_error(FATAL, "ERROR: stderr() did not return the expected value")
  end if

end program test_stderr
