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
!! Test that the call to mpp_mem_dump is functional.  On a supported OS, the return value
!! of mpp_mem_dump should be a positive integer.
program test_mpp_mem_dump

  use mpp_memutils_mod, only: mpp_mem_dump
  use platform_mod
  implicit none

  real :: memuse

  call mpp_mem_dump(memuse)

  if (memuse .LT. 0.0) then
    stop "mpp_mem_dump return negative value for memory"
  endif
end program test_mpp_mem_dump
