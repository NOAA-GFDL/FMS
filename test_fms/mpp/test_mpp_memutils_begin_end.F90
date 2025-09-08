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
!! Test that the mpp_memutils_mod module's begin and end subroutines are functional.
!! This test will exit with status zero if successful.  The script launching this test
!! program may, if desired, check if the memory output numbers are sane.
program test_mpp_memutils_init_end

  use mpp_mod, only : mpp_init, mpp_exit
  use mpp_memutils_mod, only: mpp_memuse_begin, mpp_memuse_end
  use platform_mod
  implicit none

  real, dimension(:), allocatable :: ralloc_mem

  call mpp_init()
  call mpp_memuse_begin()
  allocate(ralloc_mem(10000000))
  ralloc_mem = 0.0
  call mpp_memuse_end("Test mpp_memuse_mod")
  call mpp_exit()
end program test_mpp_memutils_init_end
