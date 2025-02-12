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
!! @brief unit tests forr diag_integral_mod
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program test_diag_integral tests mainly sum_field_2d/3d/wght_3d/hemi
!! subroutines in diag_integral_mod.  Within all the public routines in diag_integral_mod,
!! there are no APIs to retrieve the computed results.  Instead ,the results have to be
!! written out to the diag_integral.out file and these values then have to be read in
!! to check for correctness.
!! TODO:  More comprehensive tests that use more complicated data.  All the arrays are initialized
!!        to 1.0_lkind in this test.  The computed results are thus simple values with "a lot of trailing zeros"
!!        Thus, the computed values should exactly equal the expected answers.
!! TODO:  This test also only checks for one timepoint.  Additional tests that checks the evolution of the average/sum
!!        wrt time are required.


program test_diag_integral

  use diag_integral_mod
  use fms_mod, only : fms_init, fms_end
  use mpp_mod, only: mpp_init, mpp_sync, mpp_npes, mpp_get_current_pelist, mpp_error, FATAL
  use time_manager_mod, only: time_type, set_time
  use platform_mod, only : r4_kind, r8_kind

  implicit none

  integer, parameter :: nxy=20     !> number of x and y coordinate points
  integer, parameter :: nxyp=nxy+1 !> dimension of lat and lon arrays

  character(9), parameter :: field_name2='immadeup2' !> made up field name to test sum_field_2d
  character(9), parameter :: field_name3='immadeup3' !> made up field name to test sum_field_3d
  character(9), parameter :: field_namew='immadeupw' !> made up field name to test sum_field_wght_3d
  character(9), parameter :: field_nameh='immadeuph' !> made up field name to test
  character(8), parameter :: std_digits   = 'e23.15e3' !> write out precision for r8_kind data

  real(TEST_DI_KIND_) :: immadeup2(nxy,nxy)     !> array to test sum_field_2d
  real(TEST_DI_KIND_) :: immadeup3(nxy,nxy,nxy) !> array to test sum_field_3d
  real(TEST_DI_KIND_) :: immadeupw(nxy,nxy,nxy) !> array to test sum_field_wght_3d
  real(TEST_DI_KIND_) :: weight(nxy,nxy,nxy)    !> weights required to test sum_field_wght_3d
  real(TEST_DI_KIND_) :: immadeuph(nxy,nxy)     !> array to test sum_field_2d_hemi

  real(r8_kind) :: lat(nxyp,nxyp), lon(nxyp,nxyp)
  real(r8_kind) :: area(nxy,nxy)
  type(time_type) :: Time_init, Time

  !testing and generating answers
  integer :: i, j, k !> counters for do loop
  real(r8_kind) :: area_sum !> global area.  sum of the grid cell areas.
  real(r8_kind) :: itime    !> made up time
  !> The field_avg* values are only declared as r8_kind because they correspond to the values
  !! stored in field_sum, which is declared as r8_kind in diag_integral_mod
  real(r8_kind) :: field_avg2 !> result from sum_field_2d
  real(r8_kind) :: field_avg3 !> result from sum_field_3d
  real(r8_kind) :: field_avgw !> result from sum_field_avgw
  real(r8_kind) :: field_avgh !> result from field_avgh

  integer, parameter :: lkind=TEST_DI_KIND_ !> local version of TEST_DI_KIND_

  call fms_init
  call initialize_arrays

  !> The diag_integral_nml is written out by test_diag_integral2.sh

  !> Made up time.  Set the initial time.
  Time=set_time(0,1,0)

  call test_diag_integral_init       !< Test that diag_integral_init call works
  call test_diag_integral_field_init !< Register the fields
  call test_call_diag_integral_field !< call test_sum_field_* to generate and write out data to diag_integral.out
  call diag_integral_end(Time)
  call fms_end
  !> diag_integral.out file is printed out by diag_integral_end

  !> compute total area
  area_sum=real(nxy*nxy,r8_kind)

  call read_diag_integral_file      !< read in computed values in file diag_integral.out
  call test_sum_diag_integral_field !< compare read in values to the expected values.

contains
  !-------------------------------------
  !-------------------------------------
  subroutine test_diag_integral_init

    !> This test subroutine ony checks that diag_integral_init can be called successfully.

    implicit none

    call diag_integral_init(Time_init=Time_init, Time=Time, blon=lon, blat=lat, area_in=area)

  end subroutine test_diag_integral_init
  !-------------------------------------
  !-------------------------------------
  subroutine test_diag_integral_field_init

    !> This subroutine only checks that diag_integral_field_init can be called successfully.
    !> This subroutine is also required to register the fields.

    implicit none

    call diag_integral_field_init(field_name2, std_digits)
    call diag_integral_field_init(field_name3, std_digits)
    call diag_integral_field_init(field_namew, std_digits)
    call diag_integral_field_init(field_nameh, std_digits)

  end subroutine test_diag_integral_field_init
  !-------------------------------------
  !-------------------------------------
  subroutine test_call_diag_integral_field

    !> This subroutine calls sum_diag_integral_field_*
    !! The computed results are written out to diag_integral.out file (via diag_integral_output)

    implicit none
    integer :: is, ie, js, je

    is=1 ; ie=nxy
    js=1 ; je=js

    Time=set_time(0,2,0)
    call sum_diag_integral_field(field_name2, immadeup2)
    call sum_diag_integral_field(field_name3, immadeup3)
    call sum_diag_integral_field(field_namew, immadeupw)
    do js=1, nxy
       je=js
       !> sum_diag_integral_hemi sums over a slice of the sphere at a time
       call sum_diag_integral_field(field_nameh, immadeuph, is, ie, js, je)
    end do
    call diag_integral_output(Time)

  end subroutine test_call_diag_integral_field
  !-------------------------------------
  !-------------------------------------
  subroutine test_sum_diag_integral_field

    !> In diag_integral_mod, the field_averages are stored in the array field_sum
    !! as r8_kind values.  Thus, the answers here are declared as r8_kind and the
    !! comparison is between r8_kind computed results and expected answers.

    implicit none

    real(r8_kind) :: answer2, answer3, answerw, answerh

    !> expected answers
    answer2=real(nxy*nxy,r8_kind)/area_sum
    answer3=real(nxy*nxy*nxy,r8_kind)/area_sum
    answerw=answer3
    answerh=answer2/2.0_lkind

    call check_answers(answer2, field_avg2, 'sum_diag_integral_field failed for 2d')
    call check_answers(answer3, field_avg3, 'sum_diag_integral_field failed for 3d')
    call check_answers(answerw, field_avgw, 'sum_diag_integral_field failed for weighted')
    call check_answers(answerh, field_avgh, 'sum_diag_integral_field_hemi failed')

  end subroutine test_sum_diag_integral_field
  !-------------------------------------
  !-------------------------------------
  subroutine read_diag_integral_file

    character(*), parameter :: di_file='diag_integral.out'
    integer :: iunit
    character(100) :: cline1, cline2, cline3, cline4, cline5, clin6

    !> read in computed values
    open(newunit=iunit, file=di_file)
    read(iunit,*) cline1, cline2, cline3, cline4, cline5, clin6
    read(iunit,*) itime, field_avg2, field_avg3, field_avgw, field_avgh
    close(iunit)

  end subroutine read_diag_integral_file
  !-------------------------------------
  !-------------------------------------
  subroutine check_answers(answer, outresult, whoami)

    implicit none
    real(r8_kind), intent(in) :: answer, outresult
    character(len=*), intent(in) :: whoami

    if(answer.ne.outresult) then
       write(*,*) '*******************************************'
       write(*,*) 'expected', answer, 'but computed',  outresult
       call mpp_error(FATAL,'ERROR: '//trim(whoami))
    end if

  end subroutine check_answers
  !-------------------------------------
  !-------------------------------------
  subroutine initialize_arrays

    !> made up numbers

    implicit none

    lon=1.0_lkind
    lat=1.0_lkind
    area=1.0_lkind
    immadeup2=1.0_lkind
    immadeup3=1.0_lkind
    immadeupw=1.0_lkind
    immadeuph=1.0_lkind
    weight=1.0_lkind

  end subroutine initialize_arrays
  !-------------------------------------
  !-------------------------------------
end program test_diag_integral
