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

  integer :: in_unit !< Stores the returned standard input unit number
  character(len=128) :: err_msg !< Stores an arbitrary error message
  logical :: open5 !< Indicates whether unit 5 is open
  logical :: open100 !< Indicates whether unit 100 is open

  call mpp_init()

  inquire(unit=5, opened=open5)
  inquire(unit=100, opened=open100)

  in_unit = stdin()

  write(*,*) "in_unit is set to either 5 or 100."
  write(*,*) "MPP sets the status of this unit to open."
  write(*,*) "stdin() should get and return the value of in_unit."

  if (open5) then
    if (in_unit.eq.5) then
      write(*,*) "PASS: Standard fortran unit for input is open at 5,"
      write(*,*) "stdin() correctly reflects that."
    else
      err_msg = "ERROR: The unit 5 is open, but 5 is not returned by stdin()."
      call mpp_error(FATAL, err_msg)
    end if
  end if

  if (open100) then
    if (in_unit.eq.100) then
      write(*,*) "PASS: Standard fortran unit for input is open at 100,"
      write(*,*) "stdin() correctly reflects that."
    else
      err_msg = "ERROR: The unit 100 is open, but 100 is not returned by stdin()."
      call mpp_error(FATAL, err_msg)
    end if
  end if

  if (.not.(open5.or.open100)) then
    err_msg = "ERROR: Neither unit 5 nor 100 is open, &
               &expected one of these to be open."
    call mpp_error(FATAL, err_msg)
  end if

  call mpp_exit()

end program test_stdin
