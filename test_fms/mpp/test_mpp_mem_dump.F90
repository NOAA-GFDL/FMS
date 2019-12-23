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
program test_mpp_mem_dump
  use mpp_memutils_mod, only: mpp_mem_dump
#include "../../include/fms_platform.h"
  implicit none

  real :: memuse
  integer(kind=INT_KIND) :: return_val = 0

  call mpp_mem_dump(memuse)

  if (memuse .LT. 0.0) then
    return_val = 1
  endif

  stop return_val
end program test_mpp_mem_dump
