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
!! Test that the mpp_memutils_mod module's begin and end subroutines are functional.
!! This test will exit with status zero if successful.  The script launching this test
!! program may, if desired, check if the memory output numbers are sane.
program test_mpp_memutils_init_end
#include "../../include/fms_platform.h"

  use mpp_mod, only : mpp_init, mpp_exit
  use mpp_memutils_mod, only: mpp_memuse_begin, mpp_memuse_end
  implicit none

  real, dimension(:), allocatable :: ralloc_mem

  call mpp_init()
  call mpp_memuse_begin()
  allocate(ralloc_mem(10000000))
  ralloc_mem = 0.0
  call mpp_memuse_end("Test mpp_memuse_mod")
  call mpp_exit()
end program test_mpp_memutils_init_end
