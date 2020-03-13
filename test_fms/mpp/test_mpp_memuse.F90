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
!! Test case to ensure getpeakrss() is functional.  If the system is supported, the return
!! on the call to getpeakrss() will be greater than 0.
program test_mpp_memuse
#include "../../include/fms_platform.h"
  implicit none

  interface
    integer(KIND=c_size_t) function getpeakrss() result(gp) bind(c, name="getpeakrss")
      use, intrinsic :: iso_c_binding
    end function getpeakrss
  end interface

  ! Real to hold results from getpeakrss
  real :: memuse

  ! Cast to real, this is how getpeakrss is called in mpp_memutils.F90
  memuse = real(getpeakrss())

  if ( memuse < 0 ) then
    stop "getpeakrss returned negative value for memory"
  endif
end program test_mpp_memuse
