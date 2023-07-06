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
!! @description This program tests mainly the lookup* procedures in sat_vapor_pres_mod.
!!              The compute_tables, compute_es_k, and compute_es_liq_k subroutines found in
!!              this file are copied from sat_vapor_pres_k.F90 in order to generate answers.
!! TODO:  A more comprehensive testing suite for the subroutine compute_qd
!! TODO:  A more comprehensive testing suite for the compute_mrs
!! TODO:  A test to check computation of TABLE and DTABLE when do_simple=.true.
!!        (see subroutine sat_vapor_press_k_init)
!! TODO:  Testing suite to test computations involving D2TABLE, D2TABLE2, D2TABLE3
!!        Current tests for the lookup* subroutines only checks esat and desat for temperatures
!!        = TCMIN and TCMAX.  The D2* tables are not involved in the computation for these two cases.
!!        Thus, testing suite that involves these D2* table values should be incorporated here.
!! TODO:  Test the computation of nbads and test to to see if the expected error occurs if the temperature
!!        is less than TCMIN or higher than TCMAX

program test_sat_vap_pressure

use            fms_mod, only: fms_init, fms_end
use            mpp_mod, only: mpp_error, FATAL
use       platform_mod, only: r4_kind, r8_kind
use      constants_mod, only: RDGAS, RVGAS, TFREEZE
use sat_vapor_pres_mod, only: TCMIN, TCMAX, sat_vapor_pres_init, &
                              compute_qs, compute_mrs,           &
                              lookup_es,   lookup_des,  lookup_es_des,   &
                              lookup_es2,  lookup_des2, lookup_es2_des2, &
                              lookup_es3,  lookup_des3, lookup_es3_des3

implicit none

integer, parameter :: ESRES=10  !> taken from sat_vapor_pres_mod
real(r8_kind), dimension(:), allocatable :: TABLE, DTABLE, TABLE2, DTABLE2, TABLE3, DTABLE3
integer :: io, N

integer, parameter :: nml_unit_var=100
character(100) :: nml_file
logical :: test1, test2, test3, test4, test5
NAMELIST / test_sat_vapor_pres_nml/ test1, test2, test3, test4, test5

N=(TCMAX-TCMIN)*ESRES+1
allocate( TABLE(N),DTABLE(N),TABLE2(N),DTABLE2(N),TABLE3(N),DTABLE3(N) )

call fms_init()
call sat_vapor_pres_init()  !> compute tables to be used for testing
call compute_tables()       !> compute tables to generate answers/reference values

nml_file='test_sat_vapor_pres.nml'
open(unit=nml_unit_var, file=trim(nml_file), action='read')
read(unit=nml_unit_var, nml=test_sat_vapor_pres_nml,iostat=io)
close(nml_unit_var)

!CALL TESTS
if(test1) then
   write(*,*)'***TEST COMPUTE_QS 1D-3D***'
   call test_compute_qs()
end if
if(test2) then
   write(*,*)'***TEST COMPUTE_MRS 1D-3D***'
   call test_compute_mrs()
end if
if(test3) then
   write(*,*)'***TEST LOOKUP_ES,  LOOKUP_DES,  LOOKUP_ES_DES, 1D-3D***'
   call test_lookup_es_des()
end if
if(test4) then
   write(*,*)'***TEST LOOKUP_ES2, LOOKUP_DES2, LOOKUP_ES2_DES2, 1D-3D***'
   call test_lookup_es2_des2()
end if
if(test5) then
   write(*,*)'***TEST_LOOKUP_ES3, LOOKUP_DES3, LOOKUP_ES3_DES3, 1D-3D***'
   call test_lookup_es3_des3()
end if

call fms_end()

contains
  !-----------------------------------------------------------------------
  subroutine test_compute_qs()

    !> TEST:  The qsat value should equal RDGAS/RVGAS as pressure is (hypothetically) zero.
    !! The tests for this section is not comprehensive and more tests should be added.

    implicit none

    real(kind=TEST_SVP_KIND_) :: temp, press, answer, qsat
    real(kind=TEST_SVP_KIND_), dimension(1)     :: temp_1d, press_1d, answer_1d, qsat_1d
    real(kind=TEST_SVP_KIND_), dimension(1,1)   :: temp_2d, press_2d, answer_2d, qsat_2d
    real(kind=TEST_SVP_KIND_), dimension(1,1,1) :: temp_3d, press_3d, answer_3d, qsat_3d

    real(kind=r8_kind), parameter :: EPSILO=real(RDGAS,r8_kind)/real(RVGAS, r8_kind)
    integer, parameter :: lkind=TEST_SVP_KIND_ !< local kind value; using TEST_SVP_KIND_ in cases
                                               !! such as 1.0_TEST_SVP_KIND_ cannot be compiled with
                                               !! with gcc compilers.

    !---- 0d ----!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp = 270.0_lkind ; press = 0.0_lkind ; answer=real(EPSILO,lkind)
    call compute_qs(temp, press, qsat)
    call check_answer_0d( answer, qsat, 'test_compute_qs_0d')

    !---- 1d ----!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp_1d = 270.0_lkind ; press_1d = 0.0_lkind ; answer_1d=real(EPSILO,lkind)
    call compute_qs(temp_1d, press_1d, qsat_1d)
    call check_answer_1d( answer_1d, qsat_1d, 'test_compute_qs_1d')

    !---- 2d ----!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp_2d = 270.0_lkind ; press_2d = 0.0_lkind ; answer_2d=real(EPSILO,lkind)
    call compute_qs(temp_2d, press_2d, qsat_2d)
    call check_answer_2d( answer_2d, qsat_2d, 'test_compute_qs_2d')

    !---- 3d ----!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp_3d = 270.0_lkind ; press_3d = 0.0_lkind ; answer_3d=real(EPSILO,lkind)
    call compute_qs(temp_3d, press_3d, qsat_3d)
    call check_answer_3d( answer_3d, qsat_3d, 'test_compute_qs_3d')

  end subroutine test_compute_qs
  !-----------------------------------------------------------------------
  subroutine test_compute_mrs()

    !> TEST:  The qsat value should equal RDGAS/RVGAS as pressure is (hypothetically) zero.
    !! The tests for this section is not comprehensive and more tests should be added.

    implicit none
    real(kind=TEST_SVP_KIND_) :: temp, press, answer, mrsat
    real(kind=TEST_SVP_KIND_), dimension(1)     :: temp_1d, press_1d, answer_1d, mrsat_1d
    real(kind=TEST_SVP_KIND_), dimension(1,1)   :: temp_2d, press_2d, answer_2d, mrsat_2d
    real(kind=TEST_SVP_KIND_), dimension(1,1,1) :: temp_3d, press_3d, answer_3d, mrsat_3d

    real(kind=r8_kind), parameter :: EPSILO=real(RDGAS,r8_kind)/real(RVGAS, r8_kind)
    integer, parameter :: lkind=TEST_SVP_KIND_ !< local kind value; using TEST_SVP_KIND_ in cases
                                               !! such as 1.0_TEST_SVP_KIND_ cannot be compiled with
                                               !! with gcc compilers.

    !--------0d--------!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp= 270.0_lkind ; press= 0.0_lkind ; answer=real(EPSILO,lkind)
    call compute_mrs(temp, press, mrsat)
    call check_answer_0d(answer,mrsat,'test_compute_mrs_0d precision=TEST_SVP_KIND_')

    !--------1d--------!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp_1d = 270.0_lkind ; press_1d = 0.0_lkind ; answer_1d=real(EPSILO,lkind)
    call compute_mrs(temp_1d, press_1d, mrsat_1d)
    call check_answer_1d(answer_1d,mrsat_1d,'test_compute_mrs_1d precision=TEST_SVP_KIND_')

    !--------2d--------!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp_2d = 270.0_lkind ; press_2d = 0.0_lkind ; answer_2d=real(EPSILO,lkind)
    call compute_mrs(temp_2d, press_2d, mrsat_2d)
    call check_answer_2d(answer_2d,mrsat_2d,'test_compute_mrs_2d precision=TEST_SVP_KIND_')

    !--------3d--------!
    !> press is 0.  Therefore the answer should be eps=EPSILO=RDGAS/RVGAS
    temp_3d = 270.0_lkind ; press_3d = 0.0_lkind ; answer_3d=real(EPSILO,lkind)
    call compute_mrs(temp_3d, press_3d, mrsat_3d)
    call check_answer_3d(answer_3d,mrsat_3d,'test_compute_mrs_3d precision=TEST_SVP_KIND_')

  end subroutine test_compute_mrs
  !-----------------------------------------------------------------------
  subroutine test_lookup_es_des

    !> TEST:  at the minimum temperature (TCMIN), the pressures should correspond to the first element in the (D)TABLE
    !! TEST:  at the maximum temperature (TCMAX), the pressures should correspond to the last element in the (D)TABLE

    implicit none
    real(kind=TEST_SVP_KIND_) ::                   temp,    esat,    desat,    esat_answer,    desat_answer
    real(kind=TEST_SVP_KIND_), dimension(1)     :: temp_1d, esat_1d, desat_1d, esat_answer_1d, desat_answer_1d
    real(kind=TEST_SVP_KIND_), dimension(1,1)   :: temp_2d, esat_2d, desat_2d, esat_answer_2d, desat_answer_2d
    real(kind=TEST_SVP_KIND_), dimension(1,1,1) :: temp_3d, esat_3d, desat_3d, esat_answer_3d, desat_answer_3d

    integer, parameter :: lkind=TEST_SVP_KIND_ !< local kind value; using TEST_SVP_KIND_ in cases
                                               !! such as 1.0_TEST_SVP_KIND_ cannot be compiled with
                                               !! with gcc compilers

    !-----0d test-------!
    !> test lookup_es
    !! at temp=TCMIN, the answers should be TABLE(1)
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE(1), lkind)
    call lookup_es(temp,esat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es_0d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N)
    temp=real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE(N),lkind)
    call lookup_es(temp,esat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es_0d precision TCMAX')

    !> test lookup_des
    !! at temp=TCMIN, the answers should be DTABLE(1)
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer=real(DTABLE(1), lkind)
    call lookup_des(temp,desat)
    call check_answer_0d(desat_answer, desat, 'test_lookup_des_0d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE(N)
    temp=real(TCMAX,lkind)+real(TFREEZE,lkind)
    desat_answer = real(DTABLE(N),lkind)
    call lookup_des(temp,desat)
    call check_answer_0d(desat_answer, desat,   'test_lookup_es_0d TCMAX')

    !> test lookup_es_des
    !! at temp=TCMIN, the answers should be TABLE(1) and DTABLE(1) respectively
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE(1), lkind)
    desat_answer = real(DTABLE(1), lkind)
    esat = 0._lkind ; desat = 0.0_lkind
    call lookup_es_des(temp,esat,desat)
    call check_answer_0d(esat_answer,  esat,  'test_lookup_es_des_0d TCMIN')
    call check_answer_0d(desat_answer, desat, 'test_lookup_es_des_0d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N), DTABLE(N) respectively
    temp = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE(N), lkind)
    desat_answer = real(DTABLE(N), lkind)
    esat = 0._lkind ; desat = 0.0_lkind
    call lookup_es_des(temp,esat,desat)
    call check_answer_0d(esat_answer,  esat,  'test_lookup_es_des_0d TCMAX')
    call check_answer_0d(desat_answer, desat, 'test_lookup_es_des_0d TCMAX')


    !-----1d test-------!
    !> test lookup_es
    !! at temp=TCMIN, the answers should be TABLE(1)
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE(1)
    call lookup_es(temp_1d,esat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es_1d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N)
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE(N)
    call lookup_es(temp_1d,esat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es_1d TCMAX')

    !> test lookup_des
    !! at temp=TCMIN, the answers should be DTABLE(1)
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_1d = DTABLE(1)
    call lookup_des(temp_1d,desat_1d)
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_des_1d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE(N)
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_1d = DTABLE(N)
    call lookup_des(temp_1d,desat_1d)
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_des_1d TCMAX')

    !> test lookup_es_des
    !! at temp=TCMIN, the answers should be TABLE(1) and DTABLE(1) respectively
    esat_1d = 0._lkind ; desat_1d = 0._lkind
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE(1)
    desat_answer_1d = DTABLE(1)
    call lookup_es_des(temp_1d,esat_1d,desat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es_des_1d TCMIN')
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_es_des_1d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N) and DTABLE(N) respectively
    esat_1d = 0._lkind ; desat_1d = 0._lkind
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE(N)
    desat_answer_1d = DTABLE(N)
    call lookup_es_des(temp_1d,esat_1d,desat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es_des_1d TCMAX')
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_es_des_1d TCMAX')

    !-----2d test-------!
    !> test lookup_es
    !! at temp=TCMIN, the answers should be TABLE(1)
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = real(TABLE(1),lkind)
    call lookup_es(temp_2d,esat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es_2d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N)
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = real(TABLE(N),lkind)
    call lookup_es(temp_2d,esat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es_2d TCMAX')

    !> test lookup_des
    !! at temp=TCMIN, the answers should be DTABLE(1)
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_2d = DTABLE(1)
    call lookup_des(temp_2d,desat_2d)
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_des_2d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE(N)
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_2d = DTABLE(N)
    call lookup_des(temp_2d,desat_2d)
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_des_2d TCMAX')

    !> test lookup_es_des
    !! at temp=TCMIN, the answers should be TABLE(1) and DTABLE(1) respectively
    esat_2d = 0._lkind ; desat_2d = 0._lkind
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = TABLE(1)
    desat_answer_2d = DTABLE(1)
    call lookup_es_des(temp_2d,esat_2d,desat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es_des_2d TCMIN')
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_es_des_2d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N) and DTABLE(N) respectively
    esat_2d = 0._lkind ; desat_2d = 0._lkind
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = TABLE(N)
    desat_answer_2d = DTABLE(N)
    call lookup_es_des(temp_2d,esat_2d,desat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es_des_2d TCMAX')
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_es_des_2d TCMAX')

    !-----3d test-------!
    !> test lookup_es
    !! at temp=TCMIN, the answers should be TABLE(1)
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE(1)
    call lookup_es(temp_3d,esat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es_3d precision TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N)
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE(N)
    call lookup_es(temp_3d,esat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es_3d TCMAX')

    !> test lookup_des
    !! at temp=TCMIN, the answers should be DTABLE(1)
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_3d = DTABLE(1)
    call lookup_des(temp_3d,desat_3d)
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_des_3d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE(N)
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_3d = DTABLE(N)
    call lookup_des(temp_3d,desat_3d)
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_des_3d TCMAX')

    !> test lookup_es_des
    !! at temp=TCMIN, the answers should be TABLE(1) and DTABLE(1) respectively
    esat_3d = 0._lkind ; desat_3d = 0._lkind
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE(1)
    desat_answer_3d = DTABLE(1)
    call lookup_es_des(temp_3d,esat_3d,desat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es_des_3d TCMIN')
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_es_des_3d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE(N) and DTABLE(N) respectively
    esat_3d = 0._lkind ; desat_3d = 0._lkind
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE(N)
    desat_answer_3d = DTABLE(N)
    call lookup_es_des(temp_3d,esat_3d,desat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es_des_3d TCMAX')
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_es_des_3d TCMAX')

  end subroutine test_lookup_es_des
  !----------------------------------------------------------------------
  subroutine test_lookup_es2_des2

    !> TEST:  at the minimum temperature (TCMIN), the pressures should correspond to the first element in the (D)TABLE2
    !! TEST:  at the maximum temperature (TCMAX), the pressures should correspond to the last element in the (D)TABLE2

    implicit none
    real(kind=TEST_SVP_KIND_)                   :: temp,    esat,    desat,    esat_answer,    desat_answer
    real(kind=TEST_SVP_KIND_), dimension(1)     :: temp_1d, esat_1d, desat_1d, esat_answer_1d, desat_answer_1d
    real(kind=TEST_SVP_KIND_), dimension(1,1)   :: temp_2d, esat_2d, desat_2d, esat_answer_2d, desat_answer_2d
    real(kind=TEST_SVP_KIND_), dimension(1,1,1) :: temp_3d, esat_3d, desat_3d, esat_answer_3d, desat_answer_3d

    integer, parameter :: lkind=TEST_SVP_KIND_ !< local kind value; using TEST_SVP_KIND_ in cases
                                               !! such as 1.0_TEST_SVP_KIND_ cannot be compiled with
                                               !! with gcc compilers.

    !-----0d test-------!
    !> test lookup_es2
    !! at temp=TCMIN, the answers should be TABLE2(1)
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE2(1),lkind)
    call lookup_es2(temp,esat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es2_0d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N)
    temp = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE2(N),lkind)
    !! test lookup_es2
    call lookup_es2(temp,esat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es2_0d TCMAX')

    !> test lookup_des2
    !! at temp=TCMIN, the answers should be DTABLE2(1)
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer=real(DTABLE2(1),lkind)
    call lookup_des2(temp,desat)
    call check_answer_0d(desat_answer, desat, 'test_lookup_des2_0d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE2(N)
    temp = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer=real(DTABLE2(N),lkind)
    call lookup_des2(temp,desat)
    call check_answer_0d(desat_answer, desat, 'test_lookup_des2_0d TCMAX')

    !> test lookup_es2_des2
    !! at temp=TCMIN, the answers should be TABLE2(1) and DTABLE2(1) respectively
    esat = 0._lkind ; desat = 0.0_lkind
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE2(1),lkind)
    desat_answer=real(DTABLE2(1),lkind)
    call lookup_es2_des2(temp,esat,desat)
    call check_answer_0d(esat_answer,   esat, 'test_lookup_es2_des2_0d TCMIN')
    call check_answer_0d(desat_answer, desat, 'test_lookup_es2_des2_0d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N) and DTABLE2(N) respectively
    esat = 0._lkind ; desat = 0.0_lkind
    temp = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer = real(TABLE2(N),lkind)
    desat_answer=real(DTABLE2(N),lkind)
    call lookup_es2_des2(temp,esat,desat)
    call check_answer_0d(esat_answer,   esat, 'test_lookup_es2_des2_0d TCMAX')
    call check_answer_0d(desat_answer, desat, 'test_lookup_es2_des2_0d TCMAX')

    !-----1d test-------!
    !> test lookup_es2
    !! at temp=TCMIN, the answers should be TABLE2(1)
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE2(1)
    call lookup_es2(temp_1d,esat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es2_1d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N)
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE2(N)
    call lookup_es2(temp_1d,esat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es2_1d TCMAX')

    !> test lookup_des2
    !! at temp=TCMIN, the answers should be DTABLE2(1)
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_1d = DTABLE2(1)
    call lookup_des2(temp_1d,desat_1d)
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_des2_1d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE2(N)
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_1d = DTABLE2(N)
    call lookup_des2(temp_1d,desat_1d)
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_des2_1d TCMAX')

    !> test lookup_es2_des2
    !! at temp=TCMIN, the answers should be TABLE2(1) and DTABLE2(1) respectively
    esat_1d = 0._lkind ; desat_1d = 0._lkind
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE2(1)
    desat_answer_1d = DTABLE2(1)
    call lookup_es2_des2(temp_1d,esat_1d,desat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es2_des2_1d TCMIN')
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_es2_des2_1d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N) and DTABLE2(N) respectively
    esat_1d = 0._lkind ; desat_1d = 0._lkind
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = TABLE2(N)
    desat_answer_1d = DTABLE2(N)
    call lookup_es2_des2(temp_1d,esat_1d,desat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es2_des2_1d TCMAX')
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_es2_des2_1d TCMAX')


    !-----2d test-------!
    !> test lookup_es2
    !! at temp=TCMIN, the answers should be TABLE2(1)
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = TABLE2(1)
    call lookup_es2(temp_2d,esat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es2_2d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N)
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = TABLE2(N)
    call lookup_es2(temp_2d,esat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es2_2d TCMAX')

    !> test lookup_des2
    !! at temp=TCMIN, the answers should be DTABLE2(1)
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_2d = DTABLE2(1)
    call lookup_des2(temp_2d,desat_2d)
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_des2_2d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE2(N)
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_2d = DTABLE2(N)
    call lookup_des2(temp_2d,desat_2d)
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_des2_2d TCMAX')

    !> test lookup_es2_des2
    !! at temp=TCMIN, the answers should be TABLE2(1) and DTABLE2(1) respectively
    esat_2d = 0._lkind ; desat_2d = 0._lkind
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = TABLE2(1)
    desat_answer_2d = DTABLE2(1)
    call lookup_es2_des2(temp_2d,esat_2d,desat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es2_des2_2d TCMIN')
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_es2_des2_2d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N) and DTABLE2(N) respectively
    esat_2d = 0._lkind ; desat_2d = 0._lkind
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = TABLE2(N)
    desat_answer_2d = DTABLE2(N)
    call lookup_es2_des2(temp_2d,esat_2d,desat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es2_des2_2d TCMAX')
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_es2_des2_2d TCMAX')


    !-----3d test-------!
    !> test lookup_es2
    !! at temp=TCMIN, the answers should be TABLE2(1)
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE2(1)
    call lookup_es2(temp_3d,esat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es2_3d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N)
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE2(N)
    call lookup_es2(temp_3d,esat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es2_3d TCMAX')

    !> test lookup_des2
    !! at temp=TCMIN, the answers should be DTABLE2(1)
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_3d = DTABLE2(1)
    call lookup_des2(temp_3d,desat_3d)
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_des2_3d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE2(N)
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_3d = DTABLE2(N)
    call lookup_des2(temp_3d,desat_3d)
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_des2_3d TCMAX')

    !> test lookup_es2_des2
    !! at temp=TCMIN, the answers should be TABLE2(1) and DTABLE2(1) respectively
    esat_3d = 0._lkind ; desat_3d = 0._lkind
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE2(1)
    desat_answer_3d = DTABLE2(1)
    call lookup_es2_des2(temp_3d,esat_3d,desat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es2_des2_3d TCMIN')
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_es2_des2_3d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE2(N) and DTABLE2(N) respectively
    esat_3d = 0._lkind ; desat_3d = 0._lkind
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE2(N)
    desat_answer_3d = DTABLE2(N)
    call lookup_es2_des2(temp_3d,esat_3d,desat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es2_des2_3d TCMAX')
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_es2_des2_3d TCMAX')

  end subroutine test_lookup_es2_des2
  !----------------------------------------------------------------------
  subroutine test_lookup_es3_des3

    !> TEST:  at the minimum temperature (TCMIN), the pressures should correspond to the first element in the (D)TABLE3
    !! TEST:  at the maximum temperature (TCMAX), the pressures should correspond to the last element in the (D)TABLE3

    implicit none
    real(kind=TEST_SVP_KIND_)                   :: temp,    esat,    desat,    esat_answer,    desat_answer
    real(kind=TEST_SVP_KIND_), dimension(1)     :: temp_1d, esat_1d, desat_1d, esat_answer_1d, desat_answer_1d
    real(kind=TEST_SVP_KIND_), dimension(1,1)   :: temp_2d, esat_2d, desat_2d, esat_answer_2d, desat_answer_2d
    real(kind=TEST_SVP_KIND_), dimension(1,1,1) :: temp_3d, esat_3d, desat_3d, esat_answer_3d, desat_answer_3d

    integer, parameter :: lkind=TEST_SVP_KIND_ !< local kind value; using TEST_SVP_KIND_ in cases
                                               !! such as 1.0_TEST_SVP_KIND_ cannot be compiled with
                                               !! with gcc compilers.

    !-----0d test-------!
    !> test lookup_es3
    !! at temp=TCMIN, the answers should be TABLE3(1)
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer = TABLE3(1)
    call lookup_es3(temp,esat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es3_0d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N)
    temp = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer = TABLE3(N)
    call lookup_es3(temp,esat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es3_0d TCMAX')

    !> test lookup_des3
    !! at temp=TCMIN, the answers should be DTABLE3(1)
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer=DTABLE3(1)
    call lookup_des3(temp,desat)
    call check_answer_0d(desat_answer, desat, 'test_lookup_des3_0d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE3(N)
    temp = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer=DTABLE3(N)
    call lookup_des3(temp,desat)
    call check_answer_0d(desat_answer, desat, 'test_lookup_des3_0d TCMAX')

    !> test lookup_es3_des3
    !! at temp=TCMIN, the answers should be TABLE3(1) and DTABLE3(1) respectively
    esat = 0._lkind ; desat = 0.0_lkind
    temp = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer = TABLE3(1)
    desat_answer = DTABLE3(1)
    call lookup_es3_des3(temp,esat,desat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es3_des3_0d TCMIN')
    call check_answer_0d(desat_answer, desat, 'test_lookup_es3_des3_0d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N) and DTABLE3(N) respectively
    esat = 0._lkind ; desat = 0.0_lkind
    temp = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer = TABLE3(N) ; desat_answer=DTABLE3(N)
    call lookup_es3_des3(temp,esat,desat)
    call check_answer_0d(esat_answer, esat,   'test_lookup_es3_des3_0d TCMAX')
    call check_answer_0d(desat_answer, desat, 'test_lookup_es3_des3_0d TCMAX')

    !-----1d test-------!
    !> test lookup_es3
    !! at temp=TCMIN, the answers should be TABLE3(1)
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = real(TABLE3(1),lkind)
    call lookup_es3(temp_1d,esat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es3_1d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N)
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = real(TABLE3(N),lkind)
    call lookup_es3(temp_1d,esat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es3_1d TCMAX')

    !> test looup_des3
    !! at temp=TCMIN, the answers should be DTABLE3(1)
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_1d = real(DTABLE3(1),lkind)
    call lookup_des3(temp_1d,desat_1d)
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_des3_1d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE3(N)
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_1d = real(DTABLE3(N),lkind)
    call lookup_des3(temp_1d,desat_1d)
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_des3_1d TCMAX')

    !> test lookup_es3_des3
    !! at temp=TCMIN, the answers should be TABLE3(1) and DTABLE3(1) respectively
    esat_1d = 0._lkind ; desat_1d = 0._lkind
    temp_1d(1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = real(TABLE3(1),lkind)
    desat_answer_1d = real(DTABLE3(1),lkind)
    call lookup_es3_des3(temp_1d,esat_1d,desat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es3_des3_1d TCMIN')
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_es3_des3_1d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N) and DTABLE3(N) respectively
    esat_1d = 0._lkind ; desat_1d = 0._lkind
    temp_1d(1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_1d = real(TABLE3(N),lkind)
    desat_answer_1d = real(DTABLE3(N),lkind)
    call lookup_es3_des3(temp_1d,esat_1d,desat_1d)
    call check_answer_1d(esat_answer_1d, esat_1d,   'test_lookup_es3_des3_1d TCMAX')
    call check_answer_1d(desat_answer_1d, desat_1d, 'test_lookup_es3_des3_1d TCMAX')


    !-----2d test-------!
    !> test lookup_es3
    !! at temp=TCMIN, the answers should be TABLE3(1)
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = real(TABLE3(1),lkind)
    call lookup_es3(temp_2d,esat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es3_2d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N)
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = real(TABLE3(N),lkind)
    call lookup_es3(temp_2d,esat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es3_2d TCMAX')

    !> test lookup_des3
    !! at temp=TCMIN, the answers should be DTABLE3(1)
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_2d = real(DTABLE3(1),lkind)
    call lookup_des3(temp_2d,desat_2d)
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_des3_2d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE3(N)
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_2d = real(DTABLE3(N),lkind)
    call lookup_des3(temp_2d,desat_2d)
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_des3_2d TCMAX')

    !> test lookup_es3_des3
    !! at temp=TCMIN, the answers should be TABLE3(1) and DTABLE3(1) respectively
    esat_2d = 0._lkind ; desat_2d = 0._lkind
    temp_2d(1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = real(TABLE3(1),lkind)
    desat_answer_2d = real(DTABLE3(1),lkind)
    call lookup_es3_des3(temp_2d,esat_2d,desat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es3_des3_2d TCMIN')
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_es3_des3_2d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N) and DTABLE3(N) respectively
    esat_2d = 0._lkind ; desat_2d = 0._lkind
    temp_2d(1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_2d = real(TABLE3(N),lkind)
    desat_answer_2d = real(DTABLE3(N),lkind)
    call lookup_es3_des3(temp_2d,esat_2d,desat_2d)
    call check_answer_2d(esat_answer_2d, esat_2d,   'test_lookup_es3_des3_2d TCMAX')
    call check_answer_2d(desat_answer_2d, desat_2d, 'test_lookup_es3_des3_2d TCMAX')

    !-----3d test-------!
    !> test lookup_es3
    !! at temp=TCMIN, the answers should be TABLE3(1)
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE3(1)
    call lookup_es3(temp_3d,esat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es3_3d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N)
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE3(N)
    call lookup_es3(temp_3d,esat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es3_3d TCMAX')

    !> test lookup_des3
    !! at temp=TCMIN, the answers should be DTABLE3(1)
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    desat_answer_3d = DTABLE3(1)
    call lookup_des3(temp_3d,desat_3d)
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_des3_3d TCMIN')
    !! at temp=TCMAX, the answers should be DTABLE3(N)
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    desat_answer_3d = DTABLE3(N)
    call lookup_des3(temp_3d,desat_3d)
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_des3_3d TCMAX')

    !> test lookup_es3_des3
    esat_3d = 0._lkind ; desat_3d = 0._lkind
    temp_3d(1,1,1) = real(TCMIN,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE3(1)
    desat_answer_3d = DTABLE3(1)
    call lookup_es3_des3(temp_3d,esat_3d,desat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es3_des3_3d TCMIN')
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_es3_des3_3d TCMIN')
    !! at temp=TCMAX, the answers should be TABLE3(N) and DTABLE3(N) respectively
    esat_3d = 0._lkind ; desat_3d = 0._lkind
    temp_3d(1,1,1) = real(TCMAX,lkind) + real(TFREEZE,lkind)
    esat_answer_3d = TABLE3(N)
    desat_answer_3d = DTABLE3(N)
    call lookup_es3_des3(temp_3d,esat_3d,desat_3d)
    call check_answer_3d(esat_answer_3d, esat_3d,   'test_lookup_es3_des3_3d TCMAX')
    call check_answer_3d(desat_answer_3d, desat_3d, 'test_lookup_es3_des3_3d TCMAX')


  end subroutine test_lookup_es3_des3
  !----------------------------------------------------------------------
  subroutine check_answer_0d(answer,fms_result,whoami)

    implicit none
    real(TEST_SVP_KIND_), intent(in) :: answer, fms_result
    character(len=*), intent(in) :: whoami

    if(answer .ne. fms_result) then
       write(*,*) 'Expected ', answer, ' but got ', fms_result
       call mpp_error(FATAL,'ERROR:'//trim(whoami) )
    end if

  end subroutine check_answer_0d
  !-----------------------------------------------------------------------
  subroutine check_answer_1d(answer,fms_result,whoami)

    implicit none
    real(TEST_SVP_KIND_), dimension(:), intent(in) :: answer, fms_result
    character(len=*), intent(in) :: whoami

    if(answer(1) .ne. fms_result(1)) then
       write(*,*) 'Expected ', answer(1), ' but got ', fms_result(1)
       call mpp_error(FATAL,'ERROR:'//trim(whoami) )
    end if

  end subroutine check_answer_1d
  !-----------------------------------------------------------------------
  subroutine check_answer_2d(answer,fms_result,whoami)

    implicit none
    real(TEST_SVP_KIND_), dimension(:,:), intent(in) :: answer, fms_result
    character(len=*), intent(in) :: whoami

    if(answer(1,1) .ne. fms_result(1,1)) then
       write(*,*) 'Expected ', answer(1,1), ' but got ', fms_result(1,1)
       call mpp_error(FATAL,'ERROR:'//trim(whoami) )
    end if

  end subroutine check_answer_2d
  !-----------------------------------------------------------------------
  subroutine check_answer_3d(answer,fms_result,whoami)

    implicit none
    real(TEST_SVP_KIND_), dimension(:,:,:), intent(in) :: answer, fms_result
    character(len=*), intent(in) :: whoami

    if(answer(1,1,1) .ne. fms_result(1,1,1)) then
       write(*,*) 'Expected ', answer(1,1,1), ' but got ', fms_result(1,1,1)
       call mpp_error(FATAL,'ERROR:'//trim(whoami) )
    end if

  end subroutine check_answer_3d
  !-----------------------------------------------------------------------
  subroutine compute_tables

    !> This subroutine is taken from the sat_vapor_pres_init_k subroutine in sat_vapor_pres/include
    !! Thus, sat_vapor_pres_init_k subroutine is not tested and is assumed to be correct.
    !! The TABLE* and DTABLE* values are required to test compute_qs, compute_mrs, and the 3 flavors of
    !! loopup_es_des subroutines
    !! The TABLE* and DTABLE* values are computed with r8_precision.


    real(kind=r8_kind), dimension(3) :: tem, es
    real(kind=r8_kind) :: dtres, tminl, dtinvl, tepsl, tinrc, tfact
    integer :: i


    !> TCMAX, TCMIN,TFREEZE are module level variables in sat_vapor_pres_mod
    dtres  = (real(TCMAX,r8_kind)-real(TCMIN,r8_kind))/real(N-1,r8_kind)
    tminl  = real(TCMIN,r8_kind)+real(TFREEZE,r8_kind)
    dtinvl = 1.0_r8_kind/dtres
    tepsl  = 0.5_r8_kind*dtres
    tinrc  = 0.1_r8_kind*dtres
    tfact  = 5.0_r8_kind*dtinvl

    do i = 1, N
       tem(1) = tminl + dtres*real(i-1,r8_kind)
       tem(2) = tem(1)-tinrc
       tem(3) = tem(1)+tinrc
       es = compute_es_k (tem, real(TFREEZE,r8_kind))
       TABLE(i)  = es(1)
       DTABLE(i) = (es(3)-es(2))*tfact
    enddo

    do i = 1, N
        tem(1) = tminl + dtres*real(i-1,r8_kind)
        tem(2) = tem(1)-tinrc
        tem(3) = tem(1)+tinrc
        !>   pass in flag to force all values to be wrt liquid
        es = compute_es_liq_k (tem, real(TFREEZE,r8_kind))
        TABLE2(i) = es(1)
        DTABLE2(i) = (es(3)-es(2))*tfact
     enddo

    do i = 1, N
       tem(1) = tminl + dtres*real(i-1,r8_kind)
       tem(2) = tem(1)-tinrc
       tem(3) = tem(1)+tinrc
       !>   pass in flag to force all values to be wrt liquid
       es = compute_es_liq_ice_k (tem, real(TFREEZE,r8_kind))
       TABLE3(i) = es(1)
       DTABLE3(i) = (es(3)-es(2))*tfact
    enddo

  end subroutine compute_tables
  !-----------------------------------------------------------------------
  function compute_es_k(tem, TFREEZE) result (es)

    !> This subroutine is taken from the compute_es_k subroutine in sat_vapor_pres/include
    !! and is required to compute the TABLE and DTABLE  values.
    !! Thus, compute_es_k subroutine is not tested and is assumed to be correct.
    !! Since the TABLE and DTABLE values are computed with r8_precision, all variables here
    !! are in r8_kind precision.


    real(kind=r8_kind), intent(in) :: tem(:), TFREEZE
    real(kind=r8_kind) :: es(size(tem,1))

    real(kind=r8_kind)  :: x, esice, esh2o, TBASW, TBASI
    integer :: i

    real(kind=r8_kind), parameter :: ESBASW = 101324.60_r8_kind
    real(kind=r8_kind), parameter :: ESBASI = 610.71_r8_kind

    real(r8_kind), parameter :: one=1.0_r8_kind
    real(r8_kind), parameter :: ten=10.0_r8_kind

    TBASW = TFREEZE+100.0_r8_kind
    TBASI = TFREEZE
    do i = 1, size(tem)

       !>  compute es over ice

       if (tem(i) < TBASI) then
          x = -9.09718_r8_kind*(TBASI/tem(i)-one)   &
               -3.56654_r8_kind*log10(TBASI/tem(i)) &
               +0.876793_r8_kind*(one-tem(i)/TBASI) + log10(ESBASI)
          esice =ten**(x)
       else
          esice = 0.0_r8_kind
       endif

       !>  compute es over water greater than -20 c.
       !!  values over 100 c may not be valid
       !!  see smithsonian meteorological tables page 350.

       if (tem(i) > -20.0_r8_kind+TBASI) then
          x = -7.90298_r8_kind*(TBASW/tem(i)-one)   &
               +5.02808_r8_kind*log10(TBASW/tem(i)) &
               -1.3816e-07_r8_kind*(ten**((one-tem(i)/TBASW)*11.344d0)-one) &
               +8.1328e-03_r8_kind*(ten**((TBASW/tem(i)-one)*(-3.49149d0))-one) &
               +log10(ESBASW)
          esh2o = ten**(x)
       else
          esh2o = 0.0_r8_kind
       endif

       !>  derive blended es over ice and supercooled water between -20c and 0c

       if (tem(i) <= -20.0_r8_kind+TBASI) then
          es(i) = esice
       else if (tem(i) >= TBASI) then
          es(i) = esh2o
       else
          es(i) = 0.05_r8_kind*((TBASI-tem(i))*esice + (tem(i)-TBASI+20.0_r8_kind)*esh2o)
       endif

    enddo

  end function compute_es_k
!-----------------------------------------------------------------------
 function compute_es_liq_k(tem, TFREEZE) result (es)

   !> This subroutine is taken from the compute_es_liq_k subroutine in sat_vapor_pres/include
   !! and is required to compute the TABLE2 and DTABLE2 values.
   !! Thus, compute_es_liq_k subroutine is not tested and is assumed to be correct.
   !! Since the TABLE2 and DTABLE2 values are computed with r8_precision, all variables here
   !! are in r8_kind precision.

   real(kind=r8_kind), intent(in) :: tem(:), TFREEZE
   real(kind=r8_kind) :: es(size(tem,1))

   real(kind=r8_kind) :: x, esh2o, TBASW
   integer :: i

   real(kind=r8_kind), parameter :: one=1.0_r8_kind
   real(kind=r8_kind), parameter :: ten=10.0_r8_kind
   real(kind=r8_kind), parameter :: ESBASW = 101324.60_r8_kind


   TBASW = TFREEZE+100.0_r8_kind

   do i = 1, size(tem)
!>  compute es over water for all temps.
!!  values over 100 c may not be valid
!!  see smithsonian meteorological tables page 350.
      x = -7.90298_r8_kind*(TBASW/tem(i)-one) &
           +5.02808_r8_kind*log10(TBASW/tem(i)) &
           -1.3816e-07_r8_kind*(ten**((one-tem(i)/TBASW)*11.344_r8_kind)-one)    &
           +8.1328e-03_r8_kind*(ten**((TBASW/tem(i)-one)*-3.49149_r8_kind)-one)&
           +log10(ESBASW)
      esh2o = ten**(x)
      es(i) = esh2o

   enddo

 end function compute_es_liq_k
 !-----------------------------------------------------------------------
 function compute_es_liq_ice_k(tem, TFREEZE) result (es)

   !> This subroutine is taken from the compute_es_liq_ice_k subroutine in sat_vapor_pres/include
   !! and is required to compute the TABLE3 and DTABLE3 values.
   !! Thus, compute_es_liq_ice_k subroutine is not tested and is assumed to be correct.
   !! Since the TABLE3 and DTABLE3 values are computed with r8_precision, all variables here
   !! are in r8_kind precision.

   real(kind=r8_kind), intent(in) :: tem(:), TFREEZE
   real(kind=r8_kind) :: es(size(tem,1))

   real(kind=r8_kind)    :: x, TBASW, TBASI
   integer :: i

   real(kind=r8_kind), parameter :: ESBASW = 101324.60_r8_kind
   real(kind=r8_kind), parameter :: ESBASI = 610.71_r8_kind
   real(kind=r8_kind), parameter :: one=  1.0_r8_kind
   real(kind=r8_kind), parameter :: ten= 10.0_r8_kind

   TBASW = TFREEZE+100.0_r8_kind
   TBASI = TFREEZE

   do i = 1, size(tem)

     if (tem(i) < TBASI) then

!>  compute es over ice

         x = -9.09718_r8_kind*(TBASI/tem(i)-one)   &
             -3.56654_r8_kind*log10(TBASI/tem(i)) &
             +0.876793_r8_kind*(one-tem(i)/TBASI) + log10(ESBASI)
         es(i) =ten**(x)
     else

!>  compute es over water
!!  values over 100 c may not be valid
!!  see smithsonian meteorological tables page 350.

          x = -7.90298_r8_kind*(TBASW/tem(i)-one) &
              +5.02808_r8_kind*log10(TBASW/tem(i)) &
              -1.3816e-07_r8_kind*(ten**((one-tem(i)/TBASW)*11.344_r8_kind)-one)      &
              +8.1328e-03_r8_kind*(ten**((TBASW/tem(i)-one)*(-3.49149_r8_kind))-one) &
             +log10(ESBASW)
         es(i) = ten**(x)
     endif
   enddo

 end function compute_es_liq_ice_k
 !-----------------------------------------------------------------------
end program test_sat_vap_pressure
