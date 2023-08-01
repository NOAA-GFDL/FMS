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
!> @author Caitlyn McAllister
!> @brief Unit test for astronomy/daily_mean_solar interfaces
!> @email gfdl.climate.model.info@noaa.gov
!> @description Performs calculations done in astronomy_mod using the daily_mean_solar
!> interfaces using 32 and 64 bit reals
!! TODO: A more comprehensive testing suite for any daily_mean_solar_cal subroutine

program test_daily_solar

  use astronomy_mod
  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only: mpp_error, FATAL
  use time_manager_mod,  only: JULIAN, set_calendar_type, time_type, set_time
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  call fms_init
  call set_calendar_type(JULIAN)
  call astronomy_init

  call test_set_get_orbital_parameters
  call test_set_get_ref_date_of_ae
  call test_diurnal_solar
  call test_daily_mean_solar
  call test_annual_mean_solar

  call fms_end

  contains

  subroutine test_set_get_orbital_parameters

    implicit none
    real(kind=TEST_AST_KIND_) :: ecc_in, obliq_in, per_in
    real(kind=TEST_AST_KIND_) :: ecc_out, obliq_out, per_out
    integer, parameter        :: lkind = TEST_AST_KIND_

    ecc_in   = 0.0_lkind
    obliq_in = 0.0_lkind
    per_in   = 0.0_lkind

    call set_orbital_parameters(ecc_in, obliq_in, per_in)
    call get_orbital_parameters(ecc_out, obliq_out, per_out)

    call check_answers(ecc_in, ecc_out, 'test_set_get_orbital_parameters')
    call check_answers(obliq_in, obliq_out, 'test_set_get_orbital_parameters')
    call check_answers(per_in, per_out, 'test_set_get_orbital_parameters')

  end subroutine test_set_get_orbital_parameters

  !---------------------------------------------!
  subroutine test_set_get_ref_date_of_ae
  end subroutine test_set_get_ref_date_of_ae
  !---------------------------------------------!

  subroutine test_diurnal_solar

    implicit none
    real(kind=TEST_AST_KIND_), dimension(1,1) :: lat2D, lon2D, cosz2D, fracday2D
    real(kind=TEST_AST_KIND_), dimension(1)   :: lat1D, lon1D, cosz1D, fracday1D
    real(kind=TEST_AST_KIND_)                 :: lat0D, lon0D, cosz0D, fracday0D
    real(kind=TEST_AST_KIND_)                 :: gmt, time_since_ae, rrsun
    type(time_type)                           :: time_in
    integer, parameter                        :: lkind = TEST_AST_KIND_

    time_since_ae = 0.0_lkind
    gmt           = 0.0_lkind

    ! test dirunal_solar_2d
    lat2D = 0.0_lkind
    lon2D = 0.0_lkind
    call diurnal_solar(lat2D, lon2D, gmt, time_since_ae, cosz2D, fracday2D, rrsun)
    call check_answers(cosz2D(1,1), 0.0_lkind, 'test_diurnal_solar_2D cosz')
    call check_answers(fracday2D(1,1), 0.0_lkind, 'test_diurnal_solar_2D fracday')
    call check_answers(rrsun, 1.0_lkind, 'test_dirunal_solar_2D rrsun')

    ! test diurnal_solar_1d
    lat1D = 0.0_lkind
    lon1D = 0.0_lkind
    call diurnal_solar(lat1D, lon1D, gmt, time_since_ae, cosz1D, fracday1D, rrsun)
    call check_answers(cosz1D(1), 0.0_lkind, 'test_diurnal_solar_1D cosz')
    call check_answers(fracday1D(1), 0.0_lkind, 'test_diurnal_solar_1D fracday')
    call check_answers(rrsun, 1.0_lkind, 'test_dirunal_solar_1D rrsun')

    ! test diurnal_solar_0d
    lat0D = 0.0_lkind
    lon0D = 0.0_lkind
    call diurnal_solar(lat0D, lon0D, gmt, time_since_ae, cosz0D, fracday0D, rrsun)
    call check_answers(cosz0D, 0.0_lkind, 'test_diurnal_solar_0D cosz')
    call check_answers(fracday0D, 0.0_lkind, 'test_diurnal_solar_0D fracday')
    call check_answers(rrsun, 1.0_lkind, 'test_dirunal_solar_0D rrsun')

    ! test diurnal_solar_cal_2d
    lat2D = 0.0_lkind
    lon2D = 0.0_lkind
    time_in = set_time(seconds=0, days=1, ticks=0 )
    call diurnal_solar(lat2D, lon2D, time_in, cosz2D, fracday2D, rrsun)
    call diurnal_solar(lat2D, lon2D, gmt, time_since_ae, cosz2D, fracday2D, rrsun)
    call check_answers(cosz2D(1,1), 0.0_lkind, 'test_diurnal_solar_cal_2D cosz')
    call check_answers(fracday2D(1,1), 0.0_lkind, 'test_diurnal_solar_cal_2D fracday')
    call check_answers(rrsun, 1.0_lkind, 'test_dirunal_solar_cal_2D rrsun')

    ! test diurnal_solar_cal_1d
    lat1D = 0.0_lkind
    lon1D = 0.0_lkind
    time_in = set_time(seconds=0, days=1, ticks=0 )
    call diurnal_solar(lat1D, lon1D, time_in, cosz1D, fracday1D, rrsun)
    call check_answers(cosz1D(1), 0.0_lkind, 'test_diurnal_solar_cal_1D cosz')
    call check_answers(fracday1D(1), 0.0_lkind, 'test_diurnal_solar_cal_1D fracday')
    call check_answers(rrsun, 1.0_lkind, 'test_dirunal_solar_cal_1D rrsun')

    ! test diurnal_solar_cal_0d
    lat0D = 0.0_lkind
    lon0D = 0.0_lkind
    time_in = set_time(seconds=0, days=1, ticks=0 )
    call diurnal_solar(lat0D, lon0D, time_in, cosz0D, fracday0D, rrsun)
    call diurnal_solar(lat0D, lon0D, gmt, time_since_ae, cosz0D, fracday0D, rrsun)
    call check_answers(cosz0D, 0.0_lkind, 'test_diurnal_solar_cal_0D cosz')
    call check_answers(fracday0D, 0.0_lkind, 'test_diurnal_solar_cal_0D fracday')
    call check_answers(rrsun, 1.0_lkind, 'test_dirunal_solar_cal_0D rrsun')

  end subroutine test_diurnal_solar

    !---------------------------------------------!

  subroutine test_daily_mean_solar

    implicit none
    real(kind=TEST_AST_KIND_), dimension(1,1) :: lat2D, cosz2D, h_out2D
    real(kind=TEST_AST_KIND_), dimension(1)   :: lat1D, cosz1D, h_out1D, solar1D
    real(kind=TEST_AST_KIND_)                 :: lat0D, cosz0D, h_out0D
    real(kind=TEST_AST_KIND_)                 :: time_since_ae, rr_out, solar_local
    type(time_type)                           :: time_in
    integer, parameter                        :: lkind = TEST_AST_KIND_
    real(kind=TEST_AST_KIND_), parameter      :: half_pi = acos(0.0_lkind)
    real(kind=TEST_AST_KIND_), parameter :: cosz_local=1.0_lkind/half_pi
    real(kind=TEST_AST_KIND_), parameter :: hout_local=half_pi/real(PI,TEST_AST_KIND_)

    time_since_ae = 0.0_lkind
    time_in = set_time(seconds=0, days=1, ticks=0 )

    ! test daily_mean_solar_2d
    lat2D = 0.0_lkind
    call daily_mean_solar(lat2D, time_since_ae, cosz2D, h_out2D, rr_out)
    call check_answers(cosz2D(1,1), 1.0_lkind/half_pi, 'test_daily_mean_solar_2D cosz2D')
    call check_answers(h_out2D(1,1),half_pi/real(PI,TEST_AST_KIND_), 'test_diurnal_solar_2D h_out2D')
    call check_answers(rr_out, 1.0_lkind, 'test_dirunal_solar_2D rr_out')

    ! test daily_mean_solar_1d
    lat1D = 0.0_lkind
    call daily_mean_solar(lat1D, time_since_ae, cosz1D, h_out1D, rr_out)
    call check_answers(cosz1D(1), 1.0_lkind/half_pi, 'test_daily_mean_solar_1D cosz1D')
    call check_answers(h_out1D(1),half_pi/real(PI,TEST_AST_KIND_), 'test_diurnal_solar_1D h_out1D')
    call check_answers(rr_out, 1.0_lkind, 'test_dirunal_solar_1D rr_out')


    ! test daily_mean_solar_0d
    lat0D = 0.0_lkind
    call daily_mean_solar(lat0D, time_since_ae, cosz0D, h_out0D, rr_out)
    call check_answers(cosz0D, 1.0_lkind/half_pi, 'test_daily_mean_solar_0D cosz0D')
    call check_answers(h_out0D,half_pi/real(PI,TEST_AST_KIND_), 'test_diurnal_solar_0D h_out0D')
    call check_answers(rr_out, 1.0_lkind, 'test_dirunal_solar_0D rr_out')


    ! test daily_mean_solar_cal_2d
    lat2D = 0.0_lkind
    call daily_mean_solar(lat2D, time_in, cosz2D, h_out2D, rr_out)
    call check_answers(cosz2D(1,1), 1.0_lkind/half_pi, 'test_daily_mean_solar_cal_2D cosz2D')
    call check_answers(h_out2D(1,1),half_pi/real(PI,TEST_AST_KIND_), 'test_diurnal_solar_cal_2D h_out2D')
    call check_answers(rr_out, 1.0_lkind, 'test_dirunal_solar_cal_2D rr_out')

    ! test daily_mean_solar_cal_1d
    lat1D = 0.0_lkind
    call daily_mean_solar(lat1D, time_in, cosz1D, h_out1D, rr_out)
    call check_answers(cosz1D(1), 1.0_lkind/half_pi, 'test_daily_mean_solar_cal_1D cosz1D')
    call check_answers(h_out1D(1),half_pi/real(PI,TEST_AST_KIND_), 'test_diurnal_solar_cal_1D h_out1D')
    call check_answers(rr_out, 1.0_lkind, 'test_dirunal_solar_cal_1D rr_out')

    ! test daily_mean_solar_cal_0d
    lat0D = 0.0_lkind
    call daily_mean_solar(lat0D, time_in, cosz0D, h_out0D, rr_out)
    call check_answers(cosz0D, 1.0_lkind/half_pi, 'test_daily_mean_solar_cal_0D cosz0D')
    call check_answers(h_out0D,half_pi/real(PI,TEST_AST_KIND_), 'test_diurnal_solar_cal_0D h_out0D')
    call check_answers(rr_out, 1.0_lkind, 'test_dirunal_solar_cal_0D rr_out')

    ! test daily_mean_solar_2level
    lat1D = 0.0_lkind
    call daily_mean_solar(lat1D, time_since_ae, cosz1D, solar1D)
    call check_answers(cosz1D(1), cosz_local, 'test_daily_mean_solar_2level cosz2D')
    call check_answers(solar1D(1),cosz_local*hout_local, 'test_daily_mean_solar_2level solar1D')

    ! test daily_mean_solar_cal_2level
    lat1D = 0.0_lkind
    call daily_mean_solar(lat1D, time_in, cosz1D, solar1D)
    call check_answers(cosz1D(1), 1.0_lkind/half_pi, 'test_daily_mean_solar_cal_2level cosz2D')
    call check_answers(solar1D(1),cosz_local*hout_local, 'test_daily_mean_solar_cal_2level solar1D')

  end subroutine test_daily_mean_solar

  !---------------------------------------------!

  subroutine test_annual_mean_solar

    implicit none
    integer :: js, je
    real(kind=TEST_AST_KIND_), dimension(1,1) :: lat2D, solar2D, cosz2D, fracday2D
    real(kind=TEST_AST_KIND_), dimension(1)   :: lat1D, solar1D, cosz1D, fracday1D
    real(kind=TEST_AST_KIND_)                 :: rrsun
    real(kind=TEST_AST_KIND_), parameter      :: half_pi = acos(0.0_r8_kind)
    integer, parameter                        :: lkind = TEST_AST_KIND_

    js = 1 ; je = 1
    lat2D = 0.0_lkind

    call annual_mean_solar(js, je, lat2D, cosz2D, solar2D, fracday2D, rrsun)
    !if (cosz2D(1,1) .ne. 1.0_lkind/half_pi)  call mpp_error(FATAL, 'test_annual_mean_solar_2D cosz2D')
    !if (solar2D(1,1) .ne. 1.0_lkind/real(PI,TEST_AST_KIND_)) call mpp_error(FATAL, 'test_annual_mean_solar_2D solar2D')
    !if (fracday2D(1,1) .ne. 1.0_lkind/2.0_lkind)             call mpp_error(FATAL, 'test_annual_mean_solar_2D fracday2D')
    !if (rrsun .ne. 1.0_lkind)                                call mpp_error(FATAL, 'test_annual_mean_solar_2D rrsun')

    call astronomy_end
    call astronomy_init
    call test_set_get_orbital_parameters
    js = 1 ; je = 1
    lat1D = 0.0_lkind

    call annual_mean_solar(js, je, lat1D, cosz1D, solar1D, fracday1D, rrsun)
    !if (cosz1D(1) .ne. 2.0_lkind/real(PI,TEST_AST_KIND_))  call mpp_error(FATAL, 'test_annual_mean_solar_1D cosz1D')
    !if (solar1D(1) .ne. 1.0_lkind/real(PI,TEST_AST_KIND_)) call mpp_error(FATAL, 'test_annual_mean_solar_1D solar1D')
    !if (fracday1D(1) .ne. 1.0_lkind/2.0_lkind)             call mpp_error(FATAL, 'test_annual_mean_solar_1d fracday1D')
    !if (rrsun .ne. 1.0_lkind)                              call mpp_error(FATAL, 'test_annual_mean_solar_1D rrsun')

    call astronomy_end
    call astronomy_init
    call test_set_get_orbital_parameters
    lat1D = 0.0_lkind

    call annual_mean_solar(lat1D, cosz1D, solar1D)
    call test_set_get_orbital_parameters
    !if (cosz1D(1) .ne. 2.0_lkind/real(PI,TEST_AST_KIND_))  call mpp_error(FATAL, 'test_annual_mean_solar_2level cosz1D')
    !if (solar1D(1) .ne. 1.0_lkind/real(PI,TEST_AST_KIND_)) call mpp_error(FATAL, 'test_annual_mean_solar_2level solar1D')

  end subroutine test_annual_mean_solar

  !---------------------------------------------!
  subroutine check_answers( results, answers, whoami )

    implicit none
    real(TEST_AST_KIND_) :: answers, results
    character(*) :: whoami

    if (results.ne.answers) then
       write(*,*) 'EXPECTED ', answers, ' but computed ', results
       call mpp_error(FATAL, trim(whoami))
    end if


  end subroutine check_answers
  !---------------------------------------------!


end program test_daily_solar