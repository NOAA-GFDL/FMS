!> \file
!! \author @underwoo
!!
!! \section LICENSE
!!
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
!!
!! \section DESCRIPTION
!!
!! Test the mpp_print_memuse_stats is functional, printing the output to stdout.
!! For this test, the stdout() call -- which returns a Fortran unit number, is used
!! in place of a new file.
!! This test will allocate memory using two allocatable arrays, and print the memory
!! use statitics between each allocation.  The program should exit with status zero if
!! successful.  The script calling the program can check if the numbers are sane, if
!! desired.
program test_mpp_mem_print_stats_file

  use mpp_mod, only : mpp_init, mpp_exit, stdout
  use mpp_memutils_mod, only: mpp_memuse_begin, mpp_memuse_end
  use mpp_memutils_mod, only: mpp_print_memuse_stats
  use platform_mod
  implicit none

  real, dimension(:), allocatable :: ralloc_mem1, ralloc_mem2

  call mpp_init()
  call mpp_memuse_begin()
  allocate(ralloc_mem1(1000000))
  ralloc_mem1 = 0.0
  call mpp_print_memuse_stats("Test ralloc_mem1 allocated", stdout())
  allocate(ralloc_mem2(1000000))
  ralloc_mem2 = 0.0
  call mpp_print_memuse_stats("Test ralloc_mem2 allocated", stdout())
  call mpp_memuse_end("Test mpp_memuse_mod", 6)
  call mpp_exit()
end program test_mpp_mem_print_stats_file
