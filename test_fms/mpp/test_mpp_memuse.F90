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
