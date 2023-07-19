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
!! TODO: A more comprehensive testing suite with various dates and times

program test_daily_solar_cal

  use astronomy_mod
  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only: mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use time_manager_mod,  only: JULIAN, set_calendar_type, time_type, set_date
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  real(kind=r8_kind), parameter :: twopi = 2.0_r8_kind * real(PI, r8_kind)
  real(kind=r8_kind), parameter :: deg_to_rad = twopi/360.0_r8_kind

  real(kind=r8_kind)   :: ecc   = 0.01671_r8_kind !< Eccentricity of Earth's orbit [dimensionless]
  real(kind=r8_kind)   :: obliq = 23.439_r8_kind   !< Obliquity [degrees]
  real(kind=r8_kind)   :: per   = 102.932_r8_kind  !< Longitude of perihelion with respect

  type(time_type)      :: time

  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  !call test_daily_mean_solar_cal_2d
  !call test_daily_mean_solar_cal_1d
  !call test_daily_mean_solar_cal_2level
  !call test_daily_mean_solar_cal_0d

  time = set_date(2022,9,22,12,0,0)

  call fms_end()

  contains

  ! test routines below, reference values created below
  subroutine test_daily_mean_solar_cal_2d

    implicit none
    integer, parameter                        :: n = 4
    real(kind=TEST_AST_KIND_), dimension(n,n) :: lat_2d, h_2d
    real(kind=TEST_AST_KIND_)                 :: ang, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(n,n) :: cosz_2d, fracday_2d
    real(kind=TEST_AST_KIND_), dimension(n,n) :: ref_cosz_2d, ref_fracday_2d
    real(kind=TEST_AST_KIND_)                 :: rrsun, ref_rrsun, ref_dec
    integer, parameter                        :: lkind = TEST_AST_KIND_
    integer                                   :: i, j, counter

    ! set input values
    ang           = 0.0_lkind
    h_2d          = real(PI,TEST_AST_KIND_)/2.0_lkind

    counter = 1
    do i = 1, 4
      do j = 1, 4
        lat_2d(j,i) = real(counter, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
        counter = counter + 1
      end do
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, ref_dec)
    call create_ref_fracday_2d(h_2d, ref_fracday_2d)
    call create_ref_cosz_2d(lat_2d, h_2d, ref_dec, time, time_since_ae, ref_cosz_2d)

    call daily_mean_solar(lat_2d, time, cosz_2d, fracday_2d, rrsun)

    do i = 1, 4
      do j = 1, 4
        if (ref_cosz_2d(i,j) .ne. cosz_2d(i,j)) then
          print *, ref_cosz_2d(i,j), " @ ", i,j, " does not equal ", cosz_2d(i,j)
          call mpp_error(FATAL, "test_daily_mean_cal: 2d cosz value does not match reference value")
        end if

        if (ref_fracday_2d(i,j) .ne. fracday_2d(i,j)) then
          print *, ref_fracday_2d(i,j), " @ ", i,j, " does not equal ", fracday_2d(i,j)
          call mpp_error(FATAL, "test_daily_mean_cal: 2d fracday value does not match reference value")
        end if
      end do
    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " does not equal ", rrsun
      call mpp_error(FATAL, "test_daily_mean_cal: 2d rrsun value does not match reference value")
    end if

  end subroutine test_daily_mean_solar_cal_2d

  subroutine test_daily_mean_solar_cal_1d

    implicit none
    integer, parameter                      :: n = 16
    real(kind=TEST_AST_KIND_), dimension(n) :: lat_1d, h_1d
    real(kind=TEST_AST_KIND_)               :: ang, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(n) :: cosz_1d, fracday_1d
    real(kind=TEST_AST_KIND_), dimension(n) :: ref_cosz_1d, ref_fracday_1d
    real(kind=TEST_AST_KIND_)               :: rrsun, ref_rrsun, ref_dec
    integer, parameter                      :: lkind = TEST_AST_KIND_
    integer                                 :: i

    ! set input values
    ang           = 0.0_lkind
    h_1d          = real(PI,TEST_AST_KIND_)/2.0_lkind

    do i = 1, 16
      lat_1d(i) = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, ref_dec)
    call create_ref_fracday_1d(h_1d, ref_fracday_1d)
    call create_ref_cosz_1d(lat_1d, h_1d, ref_dec, time, time_since_ae, ref_cosz_1d)

    call daily_mean_solar(lat_1d, time, cosz_1d, fracday_1d, rrsun)

    do i = 1, 16
        if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
          print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
          call mpp_error(FATAL, "test_daily_mean_cal: 1d cosz value does not match reference value")
        end if

        if (ref_fracday_1d(i) .ne. fracday_1d(i)) then
          print *, ref_fracday_1d(i), " @ ", i, " does not equal ", fracday_1d(i)
          call mpp_error(FATAL, "test_daily_mean_cal: 1d fracday value does not match reference value")
        end if
    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " does not equal ", rrsun
      call mpp_error(FATAL, "test_daily_mean_cal: 1d rrsun value does not match reference value")
    end if

  end subroutine test_daily_mean_solar_cal_1d

  subroutine test_daily_mean_solar_cal_2level

    implicit none
    integer, parameter                      :: n = 16
    real(kind=TEST_AST_KIND_), dimension(n) :: lat_1d, h_1d, h_out_1d
    real(kind=TEST_AST_KIND_)               :: ang, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(n) :: cosz_1d, ref_solar, ref_fracday_1d
    real(kind=TEST_AST_KIND_), dimension(n) :: ref_cosz_1d, solar, fracday_1d
    real(kind=TEST_AST_KIND_)               :: rrsun, ref_rrsun, ref_dec
    integer, parameter                      :: lkind = TEST_AST_KIND_
    integer                                 :: i

    ! set input values
    ang           = 0.0_lkind
    h_1d          = real(PI,TEST_AST_KIND_)/2.0_lkind

    do i = 1, 16
      lat_1d(i) = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, ref_dec)
    call create_ref_fracday_1d(h_1d, ref_fracday_1d)
    call create_ref_cosz_1d(lat_1d, h_1d, ref_dec, time, time_since_ae, ref_cosz_1d)

    ! calculate ref_solar value fracday == h_out
    ref_solar = ref_cosz_1d * ref_fracday_1d * ref_rrsun

    call daily_mean_solar(lat_1d, time, cosz_1d, solar)

    do i = 1, 16
        if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
          print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
          call mpp_error(FATAL, "test_daily_mean_cal: 2level cosz value does not match reference value")
        end if

        if (ref_solar(i) .ne. solar(i)) then
          print *, ref_solar(i), " @ ", i, " does not equal ", solar(i)
          call mpp_error(FATAL, "test_daily_mean_cal: 2level solar value does not match reference value")
        end if
    end do

  end subroutine test_daily_mean_solar_cal_2level

  subroutine test_daily_mean_solar_cal_0d

    implicit none
    real(kind=TEST_AST_KIND_) :: ang, h, lat, time_since_ae
    real(kind=TEST_AST_KIND_) :: cosz, fracday, rrsun
    real(kind=TEST_AST_KIND_) :: ref_cosz, ref_rrsun, ref_dec, ref_fracday
    integer, parameter        :: lkind = TEST_AST_KIND_
    integer                   :: i

    ! set input values
    time          = set_date(2022,9,22,12,0,0)
    ang           = 0.0_lkind
    h             = real(PI,TEST_AST_KIND_)/2.0_lkind

    do i = 1, 16
      lat = real(i,TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    lat = 0.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, ref_dec)
      call create_ref_fracday_0d(h, ref_fracday)
      call create_ref_cosz_0d(lat, h, ref_dec, time, time_since_ae, ref_cosz)

      call daily_mean_solar(lat, time, cosz, fracday, rrsun)

      if (ref_cosz .ne. cosz) then
        print *, ref_cosz, "@", i, " does not equal ", cosz
        call mpp_error(FATAL, "test_daily_mean_cal: 0d cosz value does not match reference value")
      end if

      if (ref_fracday .ne. fracday) then
        print *, ref_fracday, "@", i, " does not equal ", fracday
        call mpp_error(FATAL, "test_daily_mean_cal: 0d fracday value does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, "@", i, " does not equal ", rrsun
        call mpp_error(FATAL, "test_daily_mean_cal: 0d rrsun value does not match reference value")
      end if

  end subroutine test_daily_mean_solar_cal_0d

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

  subroutine create_ref_declination(ang, dec)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: ang
    real(kind=TEST_AST_KIND_), intent(out) :: dec
    real(kind=TEST_AST_KIND_)              :: rad_obliq, sin_dec

    rad_obliq    = real(obliq, TEST_AST_KIND_) * real(deg_to_rad, TEST_AST_KIND_)
    sin_dec      = - sin(rad_obliq)*sin(ang)
    dec          = asin(sin_dec)

  end subroutine create_ref_declination

  subroutine create_ref_fracday_2d(h_2d, fracday_2d)

    implicit none
    integer, parameter                                     :: n = 4
    real(kind=TEST_AST_KIND_), intent(in), dimension(n,n)  :: h_2d
    real(kind=TEST_AST_KIND_), intent(out), dimension(n,n) :: fracday_2d
    real(kind=TEST_AST_KIND_), dimension(n,n)              :: h_out

    h_out = h_2d / real(PI, TEST_AST_KIND_)
    fracday_2d = h_out

  end subroutine create_ref_fracday_2d

  subroutine create_ref_fracday_1d(h_1d, fracday_1d)

    implicit none
    integer, parameter                                   :: n = 16
    real(kind=TEST_AST_KIND_), intent(in),  dimension(n) :: h_1d
    real(kind=TEST_AST_KIND_), intent(out), dimension(n) :: fracday_1d
    real(kind=TEST_AST_KIND_), dimension(n)              :: h_out

    h_out = h_1d / real(PI, TEST_AST_KIND_)
    fracday_1d = h_out

  end subroutine create_ref_fracday_1d

  subroutine create_ref_fracday_0d(h_0d, fracday_0d)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: h_0d
    real(kind=TEST_AST_KIND_), intent(out) :: fracday_0d
    real(kind=TEST_AST_KIND_)              :: h_out

    h_out = h_0d / real(PI, TEST_AST_KIND_)
    fracday_0d = h_out

  end subroutine create_ref_fracday_0d

  subroutine create_ref_cosz_2d(lat_in, h_in, dec_in, time_in, time_since_ae, cosz_2d)

    implicit none
    integer, parameter                                     :: n = 4
    real(kind=TEST_AST_KIND_), intent(in), dimension(n,n)  :: lat_in, h_in
    real(kind=TEST_AST_KIND_), intent(in)                  :: dec_in
    type(time_type),           intent(in)                  :: time_in
    real(kind=TEST_AST_KIND_), intent(out)                 :: time_since_ae
    real(kind=TEST_AST_KIND_), intent(out), dimension(n,n) :: cosz_2d
    integer, parameter                                     :: lkind = TEST_AST_KIND_

    ! change time type for calculation
    time_since_ae = real(orbital_time(time_in), TEST_AST_KIND_)

    ! cosz calculation
    where (h_in == 0.0_lkind)
      cosz_2d = 0.0_lkind
    else where
      cosz_2d = sin(lat_in)*sin(dec_in) + cos(lat_in)*cos(dec_in)*sin(h_in)/h_in
    end where

  end subroutine create_ref_cosz_2d

  subroutine create_ref_cosz_1d(lat_in, h_in, dec_in, time_in, time_since_ae, cosz_1d)

    implicit none
    integer, parameter                                     :: n = 16
    real(kind=TEST_AST_KIND_), intent(in), dimension(n)    :: lat_in, h_in
    real(kind=TEST_AST_KIND_), intent(in)                  :: dec_in
    type(time_type),           intent(in)                  :: time_in
    real(kind=TEST_AST_KIND_), intent(out)                 :: time_since_ae
    real(kind=TEST_AST_KIND_), intent(out), dimension(n)   :: cosz_1d
    integer, parameter                                     :: lkind = TEST_AST_KIND_

    ! change time type for calculation
    time_since_ae = real(orbital_time(time_in), TEST_AST_KIND_)

    ! cosz calculation
    where (h_in == 0.0_lkind)
      cosz_1d = 0.0_lkind
    else where
      cosz_1d = sin(lat_in)*sin(dec_in) + cos(lat_in)*cos(dec_in)*sin(h_in)/h_in
    end where

  end subroutine create_ref_cosz_1d

  subroutine create_ref_cosz_0d(lat_in, h_in, dec_in, time_in, time_since_ae, cosz_0d)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)    :: lat_in, h_in
    real(kind=TEST_AST_KIND_), intent(in)    :: dec_in
    type(time_type),           intent(in)    :: time_in
    real(kind=TEST_AST_KIND_), intent(out)   :: time_since_ae
    real(kind=TEST_AST_KIND_), intent(out)   :: cosz_0d
    integer, parameter                       :: lkind = TEST_AST_KIND_

    ! change time type for calculation
    time_since_ae = real(orbital_time(time_in), TEST_AST_KIND_)

    ! cosz calculation
    if (h_in == 0.0_lkind) then
      cosz_0d = 0.0_lkind
    else
      cosz_0d = sin(lat_in)*sin(dec_in) + cos(lat_in)*cos(dec_in)*sin(h_in)/h_in
    endif

  end subroutine create_ref_cosz_0d

end program test_daily_solar_cal