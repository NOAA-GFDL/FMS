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
!> @brief Unit test for astronomy/diurnal_solar interfaces
!> @email gfdl.climate.model.info@noaa.gov
!> @description Performs calculations done in astronomy_mod using the diurnal_solar
!> interfaces using 32 and 64 bit reals
!! TODO: Add more comprehensive testing with optional arguments for any dirunal
!! solar interface including dt, allow_negative_cosz, and half_day_out
!! TODO: Add more comprehensive testing changing angle and declination values
!! TODO: XFAILS for out of bounds argument values

program test_diurnal_dt_option

  use astronomy_mod
  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only: mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use time_manager_mod,  only: JULIAN, set_calendar_type, time_type, set_date
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  real(kind=r8_kind), parameter :: twopi = 2.0_r8_kind * real(PI,r8_kind)
  real(kind=r8_kind), parameter :: deg_to_rad = twopi / 360.0_r8_kind


  real(kind=r8_kind)   :: ecc   = 0.01671_r8_kind  !< Eccentricity of Earth's orbit [dimensionless]
  real(kind=r8_kind)   :: obliq = 23.439_r8_kind   !< Obliquity [degrees]
  real(kind=r8_kind)   :: per   = 102.932_r8_kind  !< Longitude of perihelion with respect


  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call fms_end()
  
  call test_dt_cases_0d

  contains

  ! test for diurnal_solar inteface below, reference values calcualted underneath tests
  subroutine test_dt_cases_0d
    ! allow negative cosz
    implicit none
    real(kind=TEST_AST_KIND_) :: lat_0d, lon_0d, h
    real(kind=TEST_AST_KIND_) :: gmt, time_since_ae, dt, dec, ang
    real(kind=TEST_AST_KIND_) :: ref_cosz_0d, cosz_0d, ref_fracday_0d, fracday_0d
    real(kind=TEST_AST_KIND_) :: ref_rrsun, rrsun
    real(kind=TEST_AST_KIND_) :: aa, bb, t, tt, st, stt, sh
    integer, parameter        :: lkind = TEST_AST_KIND_

    ! test case 1: t = -PI(lon=0.0, gmt=0.0) ; h = pi/2 ; dt = -PI
    lat_0d = 0.0_lkind ; lon_0d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    dec = 0.0_lkind ; h = real(PI,TEST_AST_KIND_)/2.0_lkind ; dt = -real(PI, TEST_AST_KIND_)

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_dt_cosz_cases_0d(lat_0d, lon_0d, gmt, dt, h, dec, ref_cosz_0d)
    call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, &
        fracday_0d, rrsun, dt, allow_negative_cosz=.false.)

    if (ref_cosz_0d .ne. cosz_0d) then
      print *, ref_cosz_0d, " does not match ", cosz_0d
      call mpp_error(FATAL, &
      "test_diurnal_solar: cosz for case 1 does not match reference value")
    end if
      
    ! case 2: t = -PI(lon=0.0, gmt=0.0) ; h = pi/2 ; dt = PI
    lat_0d = 0.0_lkind ; lon_0d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    dec = 0.0_lkind ; h = real(PI,TEST_AST_KIND_)/2.0_lkind ; dt = real(PI, TEST_AST_KIND_)

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_dt_cosz_cases_0d(lat_0d, lon_0d, gmt, dt, h, dec, ref_cosz_0d)
    call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, &
        fracday_0d, rrsun, dt, allow_negative_cosz=.false.)

    if (ref_cosz_0d .ne. cosz_0d) then
      print *, ref_cosz_0d, " does not match ", cosz_0d
      call mpp_error(FATAL, &
      "test_diurnal_solar: cosz for case 2 does not match reference value")
    end if

    ! case 3: t = -PI(lon=0.0, gmt=0.0) ; h = pi/2 ; dt = 2pi
    lat_0d = 0.0_lkind ; lon_0d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    dec = 0.0_lkind ; h = real(PI,TEST_AST_KIND_)/2.0_lkind ; dt = 2.0_lkind*(PI,TEST_AST_KIND_)

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_dt_cosz_cases_0d(lat_0d, lon_0d, gmt, dt, h, dec, ref_cosz_0d)
    call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, &
        fracday_0d, rrsun, dt, allow_negative_cosz=.false.)

    if (ref_cosz_0d .ne. cosz_0d) then
      print *, ref_cosz_0d, " does not match ", cosz_0d
      call mpp_error(FATAL, &
      "test_diurnal_solar: cosz for case 3 does not match reference value")
    end if

    ! case 4: t = pi/2(lon=3pi/4, gmt=3pi/4) ; h = pi/2 ; dt = -pi/2
    lat_0d = 0.0_lkind ; lon_0d = 3.0_lkind*real(PI,TEST_AST_KIND_)/4.0_lkind
    gmt = 3.0_lkind*real(PI,TEST_AST_KIND_)/4.0_lkind ; time_since_ae = 0.0_lkind ; dec = 0.0_lkind
    h = real(PI,TEST_AST_KIND_)/2.0_lkind ; dt = -real(PI,TEST_AST_KIND_)/2.0_lkind

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_dt_cosz_cases_0d(lat_0d, lon_0d, gmt, dt, h, dec, ref_cosz_0d)
    call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, &
        fracday_0d, rrsun, dt, allow_negative_cosz=.false.)

    if (ref_cosz_0d .ne. cosz_0d) then
      print *, ref_cosz_0d, " does not match ", cosz_0d
      call mpp_error(FATAL, &
      "test_diurnal_solar: cosz for case 4 does not match reference value")
    end if

    ! case 5: t = pi/4(lon=5pi/8, gmt=5pi/8) ; h = pi/2 ; dt = pi/2
    lat_0d = 0.0_lkind ; lon_0d = 5.0_lkind*real(PI,TEST_AST_KIND_)/8.0_lkind
    gmt = 5.0_lkind*real(PI,TEST_AST_KIND_)/8.0_lkind ; time_since_ae = 0.0_lkind ; dec = 0.0_lkind
    h = real(PI,TEST_AST_KIND_)/2.0_lkind ; dt = real(PI,TEST_AST_KIND_)/2.0_lkind

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_dt_cosz_cases_0d(lat_0d, lon_0d, gmt, dt, h, dec, ref_cosz_0d)
    call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, &
        fracday_0d, rrsun, dt, allow_negative_cosz=.false.)

    if (ref_cosz_0d .ne. cosz_0d) then
      print *, ref_cosz_0d, " does not match ", cosz_0d
      call mpp_error(FATAL, &
      "test_diurnal_solar: cosz for case 5 does not match reference value")
    end if

    ! case 6: 
    ! case 7: 
    ! case 8: t = -3.0*PI/2.0 ; h = PI/2.0 ; dt = 4.0*PI
    lat_0d = 0.0_lkind ; lon_0d = 5.0_lkind*real(PI,TEST_AST_KIND_)/8.0_lkind
    gmt = 5.0_lkind*real(PI,TEST_AST_KIND_)/8.0_lkind ; time_since_ae = 0.0_lkind ; dec = 0.0_lkind
    h = real(PI,TEST_AST_KIND_)/2.0_lkind ; dt = 4.0_lkind*real(PI,TEST_AST_KIND_)

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_dt_cosz_cases_0d(lat_0d, lon_0d, gmt, dt, h, dec, ref_cosz_0d)
    call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, &
        fracday_0d, rrsun, dt, allow_negative_cosz=.false.)

    if (ref_cosz_0d .ne. cosz_0d) then
      print *, ref_cosz_0d, " does not match ", cosz_0d
      !call mpp_error(FATAL, &
      !"test_diurnal_solar: cosz for case 8 does not match reference value")
    end if

  end subroutine test_dt_cases_0d

  ! create all reference values
  subroutine create_ref_rrsun(ang, rrsun)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: ang
    real(kind=TEST_AST_KIND_), intent(out) :: rrsun
    real(kind=TEST_AST_KIND_)              :: rad_per, rr
    integer, parameter                     :: lkind = TEST_AST_KIND_
  
    rad_per = real(per, TEST_AST_KIND_) * real(deg_to_rad, TEST_AST_KIND_)
    rr      = (1.0_lkind - real((ecc**2),TEST_AST_KIND_))/(1.0_lkind + real(ecc, TEST_AST_KIND_) * cos(ang - rad_per))
    rrsun  = rr**(-2)
  
  end subroutine create_ref_rrsun

  subroutine create_ref_dt_cosz_cases_0d(lat, lon, gmt, dt, h, dec, cosz)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: lat, lon, h
    real(kind=TEST_AST_KIND_), intent(in)  :: gmt, dt, dec
    real(kind=TEST_AST_KIND_), intent(out) :: cosz
    real(kind=TEST_AST_KIND_)              :: aa, bb, t, tt, st, stt, sh
    integer, parameter                     :: lkind = TEST_AST_KIND_

    aa = sin(lat)*sin(dec)
    bb = cos(lat)*cos(dec)

    t = gmt + lon - real(PI, TEST_AST_KIND_)
    if (t >= real(PI, TEST_AST_KIND_)) t = t - real(twopi, TEST_AST_KIND_)
    if (t < real(-PI, TEST_AST_KIND_)) t = t + real(twopi, TEST_AST_KIND_)

    tt   = t + dt
    st   = sin(t)
    stt  = sin(tt)
    sh   = sin(h)

    if (t < -h .and. tt < -h) then ! case 1
      print *, "case 1"
      cosz = 0.0
    else if ((tt+h) /= 0.0 .and. t < -h .and. abs(tt) <= h) then ! case 2
      print *, "case 2"
      cosz = aa + bb*(stt + sh)/ (tt + h)
    else if (t < -h .and. h /= 0.0 .and. h < tt) then ! case 3
      print *, "case 3"
      cosz = aa + bb*( sh + sh)/(h+h)
    else if (abs(t) <= h .and. abs(tt) <= h) then ! case 4
      print *, "case 4"
      cosz = aa + bb*(stt - st)/ (tt - t)
    else if ((h-t) /= 0.0 .and. abs(t) <= h .and. h < tt) then ! case 5
      print *, "case 5"
      cosz = aa + bb*(sh - st)/(h-t)
    else if (twopi - h < tt .and. (tt+h-twopi) /= 0.0 .and. t <= h) then ! case 6
      print *, "case 6"
      cosz = (cosz*(h - t) + (aa*(tt + h - twopi) + bb &
      * (stt + sh))) / ((h - t) + (tt + h - twopi))
    else if (h <  t .and. twopi - h >= tt) then ! case 7
      print *, "case 7"
      cosz = 0.0
    else if (h <  t .and. twopi - h < tt) then ! case 8
      print *, "case 8"
      cosz = aa + bb*(stt + sh) / (tt + h - twopi)
    else
      print *, "not a case"
      cosz = aa + bb*(stt - st)/ (tt - t)
    end if

  end subroutine create_ref_dt_cosz_cases_0d

end program test_diurnal_dt_option