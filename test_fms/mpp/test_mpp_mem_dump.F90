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
!! Test that the call to mpp_mem_dump is functional.  On a supported OS, the return value
!! of mpp_mem_dump should be a positive integer.
program test_mpp_mem_dump
  use mpp_memutils_mod, only: mpp_mem_dump
#include "../../include/fms_platform.h"
  implicit none

  real :: memuse

  call mpp_mem_dump(memuse)

  if (memuse .LT. 0.0) then
    stop "mpp_mem_dump return negative value for memory"
  endif
end program test_mpp_mem_dump
