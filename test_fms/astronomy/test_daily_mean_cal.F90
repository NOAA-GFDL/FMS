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

program test_daily_solar_cal

  use astronomy_mod
  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only: mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use time_manager_mod,  only: JULIAN, set_calendar_type, time_type, set_date
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  integer, parameter                   :: lkind = TEST_AST_KIND_
  real(kind=TEST_AST_KIND_), parameter :: twopi = 2.0_lkind * real(PI, TEST_AST_KIND_)
  real(kind=TEST_AST_KIND_), parameter :: deg_to_rad = twopi/360.0_lkind
  

  real(kind=TEST_AST_KIND_)   :: ecc   = 0.01671_lkind  !< Eccentricity of Earth's orbit [dimensionless]
  real(kind=TEST_AST_KIND_)   :: obliq = 23.439_lkind   !< Obliquity [degrees]
  real(kind=TEST_AST_KIND_)   :: per   = 102.932_lkind  !< Longitude of perihelion with respect

  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call test_daily_mean_solar_cal_2d
  !call test_daily_mean_solar_cal_1d
  !call test_daily_mean_cal_2level
  call test_daily_mean_solar_cal_0d

  call fms_end()

  contains

  ! create reference values to test daily_mean_solar_2d
  subroutine create_ref_values_cal_2d(lat, dec, ang, time_since_ae, h, cosz, rr_out, fracday)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in), dimension(4,4)  :: lat, h
    real(kind=TEST_AST_KIND_), intent(in)                  :: ang
    type(time_type)                                        :: time
    real(kind=TEST_AST_KIND_), intent(out)                 :: time_since_ae
    real(kind=TEST_AST_KIND_), intent(inout)               :: dec
    real(kind=TEST_AST_KIND_), intent(out), dimension(4,4) :: cosz, fracday
    real(kind=TEST_AST_KIND_), dimension(4,4)              :: h_out
    real(kind=TEST_AST_KIND_), intent(out)                 :: rr_out
    real(kind=TEST_AST_KIND_)                              :: rad_per, rr, rad_obliq, sin_dec
    integer, parameter                                     :: lkind = TEST_AST_KIND_


    ! change time type for reference time
    time_since_ae = orbital_time(time)

    ! rr_out calculation
    rad_per = per * deg_to_rad
    rr = (1.0_lkind - ecc**2)/(1.0_lkind + ecc * cos(ang - rad_per)) ! module vars need to be TEST_KIND
    rr_out = rr**(-2)

    ! declination calculation
    rad_obliq    = obliq * deg_to_rad
    sin_dec      = - sin(rad_obliq)*sin(ang)
    dec          = asin(sin_dec)
    
    ! cosz calculation
    where (h == 0.0_lkind)
          cosz = 0.0_lkind
    elsewhere
          cosz = sin(lat)*sin(dec) + cos(lat)*cos(dec)*sin(h)/h
    end where

    ! h_out
    h_out = h / real(PI, TEST_AST_KIND_)
    fracday = h_out

  end subroutine create_ref_values_cal_2d

  subroutine test_daily_mean_solar_cal_2d

    implicit none
    real(kind=TEST_AST_KIND_), dimension(4,4) :: lat
    type(time_type)                           :: time
    real(kind=TEST_AST_KIND_)                 :: rr_out, ref_rr_out
    real(kind=TEST_AST_KIND_)                 :: ang, dec, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(4,4) :: cosz, fracday, h
    real(kind=TEST_AST_KIND_), dimension(4,4) :: ref_cosz, ref_fracday
    integer                                   :: i, j, counter
    integer, parameter                        :: lkind = TEST_AST_KIND_

    time          = set_date(2022,9,22,12,0,0)         !time of year (NH autumnal equinox @ noon)
    ang           = 0.0_lkind                          !angle(time_since_ae)
    h             = real(PI,TEST_AST_KIND_)/2.0_lkind  !half_day(lat, dec)
    lat           = 0.0_lkind

    counter = 1
    do i = 1, 4
      do j = 1, 4
        lat(j,i) = real(counter, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
        counter = counter + 1
      end do
    end do

    call create_ref_values_cal_2d(lat, dec, ang, time_since_ae, h, ref_cosz, ref_rr_out, ref_fracday)
    call daily_mean_solar(lat, time, cosz, fracday, rr_out)
    print *, "2d", ref_cosz(1,1), cosz(1,1), abs(ref_cosz(1,1) - cosz(1,1))
    print *, "2d", ref_fracday(1,1), fracday(1,1), abs(ref_fracday(1,1) - fracday(1,1))
    print *, "2d", ref_rr_out, rr_out, abs(ref_rr_out - rr_out)

    !do i = 1, 4
    !  do j = 1, 4
    !    if (abs(ref_cosz(i,j) - cosz(i,J)) .gt. real(1E-07,TEST_AST_KIND_)) then
    !      print *, ref_cosz, " @ ", i,j, " does not equal ", cosz
    !      call mpp_error(FATAL, "test_daily_mean_solar: 2d cosz value does not match reference value")
    !    end if
    !    if (abs(ref_h_out(i,j) - h_out(i,j)) .gt. real(1E-07,TEST_AST_KIND_)) then
    !      print *, ref_h_out, " @ ", i,j, " does not equal ", h_out
    !      call mpp_error(FATAL, "test_daily_mean_solar: 2d h_out value does not match reference value")
    !    end if
    !  end do
    !end do

    !if (abs(ref_rr_out - rr_out) .gt. real(1E-07,TEST_AST_KIND_)) then
    !  print *, ref_rr_out, " does not equal ", rr_out
    !  call mpp_error(FATAL, "test_daily_mean_solar: 2d rr_out value does not match reference value")
    !end if

  end subroutine test_daily_mean_solar_cal_2d 

  ! create reference values to test daily_mean_solar_0d
  subroutine create_ref_values_cal_0d(lat, dec, ang, time_since_ae, h, cosz, rr_out, fracday)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: lat, h
    real(kind=TEST_AST_KIND_), intent(in)                  :: ang
    type(time_type)                                        :: time
    real(kind=TEST_AST_KIND_), intent(out)                 :: time_since_ae
    real(kind=TEST_AST_KIND_), intent(inout)               :: dec
    real(kind=TEST_AST_KIND_), intent(out) :: cosz, fracday
    real(kind=TEST_AST_KIND_)             :: h_out
    real(kind=TEST_AST_KIND_), intent(out)                 :: rr_out
    real(kind=TEST_AST_KIND_)                              :: rad_per, rr, rad_obliq, sin_dec
    integer, parameter                                     :: lkind = TEST_AST_KIND_

    ! change time type for calculation
    time_since_ae = orbital_time(time)

    ! rr_out calculation
    rad_per = per * deg_to_rad
    rr      = (1.0_lkind - ecc**2)/(1.0_lkind + ecc * cos(ang - rad_per))
    rr_out  = rr**(-2)

    ! declination calculation
    rad_obliq    = obliq * deg_to_rad
    sin_dec      = - sin(rad_obliq)*sin(ang)
    dec          = asin(sin_dec)
    
    ! cosz calculation
    if (h == 0.0_lkind) then
      cosz = 0.0_lkind
    else
      cosz = sin(lat)*sin(dec) + cos(lat)*cos(dec)*sin(h)/h
    endif

    ! h_out
    h_out = h / real(PI, TEST_AST_KIND_)
    fracday = h_out

  end subroutine create_ref_values_cal_0d
  
  subroutine test_daily_mean_solar_cal_0d

    implicit none
    real(kind=TEST_AST_KIND_) :: time_since_ae, ang, dec, h, lat
    real(kind=TEST_AST_KIND_) :: cosz, fracday, rr_out
    real(kind=TEST_AST_KIND_) :: ref_cosz, ref_rr_out, ref_fracday
    type(time_type)           :: time
    integer                   :: i
    integer, parameter        :: lkind = TEST_AST_KIND_

    
    time          = set_date(2022,9,22,12,0,0)
    ang           = 0.0_lkind
    h             = real(PI,TEST_AST_KIND_)/2.0_lkind
    lat           = 0.0_lkind

    !do i = 0, 16
    !  lat = real(i,TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call daily_mean_solar(lat, time, cosz, fracday, rr_out)
      call create_ref_values_cal_0d(lat, dec, ang, time_since_ae, h, ref_cosz, ref_rr_out, ref_fracday)

      print *, "0d", ref_cosz, cosz, abs(ref_cosz-cosz)
      print *, "0d", ref_rr_out, rr_out, abs(ref_rr_out-rr_out)
      print *, "0d", ref_fracday, fracday, abs(ref_fracday-fracday)

      !if (abs(ref_cosz - cosz) .gt. real(1E-03,TEST_AST_KIND_)) then
      !  print *, ref_cosz, "@", i, " does not equal ", cosz
      !  call mpp_error(FATAL, "test_daily_mean_solar_cal: 0d cosz value does not match reference value")
      !end if

      !if (abs(ref_rr_out - rr_out) .gt. real(1E-03,TEST_AST_KIND_)) then
      !  print *, ref_rr_out, "@", i, " does not equal ", rr_out
      !  call mpp_error(FATAL, "test_daily_mean_solar_cal: 0d rr_out value does not match reference value")
      !end if

      !if (abs(ref_fracday - fracday) .gt. real(1E-03,TEST_AST_KIND_)) then
      !  print *, ref_fracday, "@", i, " does not equal ", fracday
      !  call mpp_error(FATAL, "test_daily_mean_solar_cal: 0d fracday value does not match reference value")
      !end if
    !end do

  end subroutine test_daily_mean_solar_cal_0d

end program test_daily_solar_cal