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
!> @author Lauren Chilutti
!> @description This test program is for testing the system_clock
!! routine.  The test initializes mpp, then calls the function 
!! SYSTEM_CLOCK.  SYSTEM_CLOCK returns COUNT, COUNT_RATE, and/or 
!! COUNT_MAX (all optional). 
program test_system_clock
#include <system_clock.h>
  use mpp_mod, only : mpp_init, mpp_exit, stderr, stdout, mpp_error, FATAL
  use mpp_mod, only : input_nml_file
  implicit none

  integer :: count, count_rate, count_max, mpicount, mpicount_rate, mpicount_max

  call SYSTEM_CLOCK(count, count_rate, count_max)
  if (count < count_max) then
    call mpp_error(FATAL, "SYSTEM_CLOCK returns a count that is greater than count_max"
  endif
  write(*,*) count, count_rate, count_max

  call mpp_init()  
  call SYSTEM_CLOCK(mpicount, mpicount_rate,mpicount_max)
  write(*,*) mpicount, mpicount_rate, mpicount_max
  call mpp_exit()
end program test_system_clock
