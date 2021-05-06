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
!! @author Lauren Chilutti
!! @brief Test program for the system_clock routine.
!! @email gfdl.climate.model.info@noaa.gov
!! @description This test program is for testing the system_clock
!! routine.  The test initializes mpp, then calls the function
!! SYSTEM_CLOCK twice.  SYSTEM_CLOCK returns COUNT, COUNT_RATE, and/or
!! COUNT_MAX (all optional).  This test returns errors if COUNT,
!! COUNT_RATE, and/or COUNT_MAX are undefined, if the COUNT is
!! greater than COUNT_MAX, that the second call returns the
!! same value for COUNT_MAX and COUNT_RATE as the first call, and that
!! the second call returns a COUNT value greater than that of the
!! first call.
module include_files_mod
  use platform_mod
  logical :: first_call_system_clock_mpi=.TRUE.
contains
#include "../../mpp/include/system_clock.h"
end module include_files_mod

program test_system_clock
  use include_files_mod
  use mpp_mod, only : mpp_init, mpp_init_test_init_true_only, stderr, stdout, mpp_error, FATAL
  implicit none

  integer(i8_kind) :: count1, count_rate1, count_max1, count2, count_rate2, count_max2
  integer :: ierr
!> Initialize mpp
  call mpp_init(test_level=mpp_init_test_init_true_only)
!> Call system_clock and ensure output is not undefined
  call system_clock(count1, count_rate1, count_max1)
  write(*,*) count1, count_rate1, count_max1
  call system_clock(count2, count_rate2, count_max2)
  write(*,*) count2, count_rate2, count_max2
!> Check that the count is less than the count_max
  if (count1 .gt. count_max1) then
    call mpp_error(FATAL, "SYSTEM_CLOCK returns a count that is greater than count_max")
  endif
  if (count2 > count_max2) then
    call mpp_error(FATAL, "SYSTEM_CLOCK returns a count that is greater than count_max")
  endif
!> Check that the second call to system_clock function returns the same COUNT_MAX and COUNT_RATE, but a larger COUNT
  if (count1 .gt. count2) then
    call mpp_error(FATAL, "count2 is less than count1")
  endif
  if (count_rate1 .ne. count_rate2) then
    call mpp_error(FATAL, "count rates are not equal")
  endif
  if (count_max1 .ne. count_max2) then
    call mpp_error(FATAL, "count maxes are not equal")
  endif
  call MPI_FINALIZE(ierr)
end program test_system_clock
