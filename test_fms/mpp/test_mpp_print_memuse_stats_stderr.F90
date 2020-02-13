!> \file
!! \author @underwoo
!!
!! \section LICENSE
!!
!***********************************************************************
!!                   GNU Lesser General Public License
!!
!! This file is part of the GFDL Flexible Modeling System (FMS).
!!
!! FMS is free software: you can redistribute it and/or modify it under
!! the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or (at
!! your option) any later version.
!!
!! FMS is distributed in the hope that it will be useful, but WITHOUT
!! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!! for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!!
!! \section DESCRIPTION
!!
!! Test the mpp_print_memuse_stats is functional, and prints to stderr, the default.
!! This test will allocate memory using two allocatable arrays, and print the memory
!! use statitics between each allocation.  The program should exit with status zero if
!! successful.  The script calling the program can check if the numbers are sane, if
!! desired.
program test_mpp_mem_print_stats_stderr
#include "../../include/fms_platform.h"

  use mpp_mod, only : mpp_init, mpp_exit
  use mpp_memutils_mod, only: mpp_memuse_begin, mpp_memuse_end
  use mpp_memutils_mod, only: mpp_print_memuse_stats
  implicit none

  real, dimension(:), allocatable :: ralloc_mem1, ralloc_mem2

  call mpp_init()
  call mpp_memuse_begin()
  allocate(ralloc_mem1(1000000))
  ralloc_mem1 = 0.0
  call mpp_print_memuse_stats("Test ralloc_mem1 allocated")
  allocate(ralloc_mem2(1000000))
  ralloc_mem2 = 0.0
  call mpp_print_memuse_stats("Test ralloc_mem2 allocated")
  call mpp_memuse_end("Test mpp_memuse_mod")
  call mpp_exit()
end program test_mpp_mem_print_stats_stderr
