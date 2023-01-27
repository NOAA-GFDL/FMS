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
!> @file
!! @brief unit test for the mpp_root_pe() function
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests some of the procedures in the sat_vap_pressure_mod module.

program test_sat_vap_pressure

use            fms_mod, only: fms_init, fms_end
use            mpp_mod, only: mpp_error, FATAL
use sat_vapor_pres_mod, only: sat_vapor_pres_init, compute_qs, compute_mrs
use       platform_mod, only: r4_kind, r8_kind
use      constants_mod, only: RDGAS, RVGAS

implicit none

call fms_init()


call sat_vapor_pres_init()
call test_compute_qs_k_0d()
call test_mrs_k_0d()

call fms_end()

contains
  !-----------------------------------------------------------------------
  subroutine test_compute_qs_k_0d()

    implicit none

    real(kind=r4_kind) :: temp4, press4, answer4, qsat4
    real(kind=r8_kind) :: temp8, press8, answer8, qsat8
    real(kind=r8_kind), parameter :: EPSILO=real(RDGAS,r8_kind)/real(RVGAS, r8_kind)

    !test 1:   press is 0.  Therefore answer should be eps=EPSILO=RDGAS/RVGAS
    temp4 = 270.0_r4_kind ; press4 = 0.0_r4_kind ; answer4=real(EPSILO,r4_kind)
    temp8 = 270.0_r8_kind ; press8 = 0.0_r8_kind ; answer8=EPSILO
    ! test r4
    call compute_qs(temp4, press4, qsat4)
    if( qsat4 .ne. answer4 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r4")
    ! test r8
    call compute_qs(temp8, press8, qsat8)
    if( qsat8 .ne. answer8 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r8")


  end subroutine test_compute_qs_k_0d
  !-----------------------------------------------------------------------
  subroutine test_mrs_k_0d()

    implicit none
    real(kind=r4_kind) :: temp4, press4, answer4, mrsat4
    real(kind=r8_kind) :: temp8, press8, answer8, mrsat8
    real(kind=r8_kind), parameter :: EPSILO=real(RDGAS,r8_kind)/real(RVGAS, r8_kind)

    !press is 0.  Therefore answer should be eps=EPSILO=RDGAS/RVGAS
    temp4 = 270.0_r4_kind ; press4 = 0.0_r4_kind ; answer4=real(EPSILO,r4_kind)
    temp8 = 270.0_r8_kind ; press8 = 0.0_r8_kind ; answer8=EPSILO
    ! test r4
    call compute_mrs(temp4, press4, mrsat4)
    if( mrsat4 .ne. answer4 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r4")
    ! test r8
    call compute_mrs(temp8, press8, mrsat8)
    if( mrsat8 .ne. answer8 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r8")

  end subroutine test_mrs_k_0d
  !-----------------------------------------------------------------------


end program test_sat_vap_pressure
