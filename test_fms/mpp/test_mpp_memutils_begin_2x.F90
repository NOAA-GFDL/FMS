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
!! This test program calls `mpp_memuse_begin` multiple times.  In the MPP
!! code this is an error, which should be caught.  This program should exit
!! non-zero.
program test_mpp_memutils_init_end

  use mpp_mod, only : mpp_init, mpp_exit
  use mpp_memutils_mod, only: mpp_memuse_begin, mpp_memuse_end
  use platform_mod
  implicit none

  real, dimension(:), allocatable :: ralloc_mem

  call mpp_init()
  ! Calling mpp_memuse_begin multiple times
  call mpp_memuse_begin()
  call mpp_memuse_begin()
  call mpp_exit()
end program test_mpp_memutils_init_end
