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
!! solar interface
!! TODO: XFAILS for out of bounds argument values
!! TODO: A more comprehensive testing suite for any dirunal_solar_cal routines

program test_diurnal_solar_cal

  use astronomy_mod
  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only: mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use time_manager_mod,  only: JULIAN, set_calendar_type, time_type, set_date
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  real(kind=TEST_AST_KIND_), parameter :: twopi = real(2.0,TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)
  real(kind=TEST_AST_KIND_), parameter :: deg_to_rad = twopi / real(360.0,TEST_AST_KIND_)
  

  real(kind=TEST_AST_KIND_)   :: ecc   = 0.01671  !< Eccentricity of Earth's orbit [dimensionless]
  real(kind=TEST_AST_KIND_)   :: obliq = 23.439   !< Obliquity [degrees]
  real(kind=TEST_AST_KIND_)   :: per   = 102.932  !< Longitude of perihelion with respect


  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call fms_end()

  call test_diurnal_solar_2d
  call test_diurnal_solar_1d
  call test_diurnal_solar_0d

  contains

  subroutine create_ref_values_0d(lat, lon, gmt, dec, ang, h, cosz, rr_out, fracday)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)    :: lat, lon, ang, h
    real(kind=TEST_AST_KIND_), intent(inout) :: dec, gmt
    real(kind=TEST_AST_KIND_), intent(out)   :: cosz, rr_out, fracday
    real(kind=TEST_AST_KIND_)                :: rad_per, rr, rad_obliq, sin_dec
    real(kind=TEST_AST_KIND_)                :: aa, bb, t
    integer, parameter                       :: lkind = TEST_AST_KIND_

    ! rr_out calculation
    rad_per = per * deg_to_rad
    rr      = (1.0_lkind - ecc**2)/(1.0_lkind + ecc * cos(ang - rad_per))
    rr_out  = rr**(-2)

    ! declination calculation
    rad_obliq    = obliq * deg_to_rad
    sin_dec      = - sin(rad_obliq)*sin(ang)
    dec          = asin(sin_dec)

    ! cosz calculation, cosz CANNOT be negative
    t = gmt + lon - real(PI, TEST_AST_KIND_)
    if (t >= real(PI, TEST_AST_KIND_))  t = t - real(twopi, TEST_AST_KIND_)
    if (t < real(-PI, TEST_AST_KIND_))  t = t + real(twopi, TEST_AST_KIND_)

    aa = sin(lat)*sin(dec)
    bb = cos(lat)*cos(dec)
    
    if (.not. Lallow_negative) then
      where (abs(t) < h)
        cosz    = aa + bb*cos(t)
        fracday = 1.0_lkind
      elsewhere
        cosz    = 0.0_lkind
        fracday = 0.0_lkind
  end where
  else
      cosz = aa + bb*cos(t)
      where (abs(t) < h)
        fracday = 1.0_lkind
      elsewhere
        fracday = 0.0_lkind
      end where
  end if

  end subroutine create_ref_values_0d

  subroutine test_diurnal_solar_0d

    implicit none
    real(kind=TEST_AST_KIND_) :: lat, lon, gmt, time_since_ae
    real(kind=TEST_AST_KIND_) :: cosz, fracday, rrsun
    real(kind=TEST_AST_KIND_) :: ref_cosz, ref_fracday, ref_rrsun
    real(kind=TEST_AST_KIND_) :: h, ang, dec
    integer                   :: i, j, k
    integer, parameter        :: lkind = TEST_AST_KIND_

    ! test by only changing lat, allow negative cosz
    lon = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI,TEST_AST_KIND_)/2.0_lkind

    lat = 0.0_lkind
    do i = 0, 16
      lat = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
      call diurnal_solar(lat, lon, gmt, time_since_ae, cosz, fracday, rrsun, allow_negative_cosz=.true.)
      call create_ref_values_0d(lat, lon, gmt, dec, ang, h, ref_cosz, ref_rrsun, ref_fracday)

      if (abs(ref_cosz - cosz) .gt. real(1E-07,TEST_AST_KIND_)) then
        print *, ref_cosz, " @ ", i, " does not equal ", cosz 
        call mpp_error(FATAL, "test_diurnal_solar: 0d cosz value does not match reference value")
      end if

      if (abs(ref_fracday - fracday) .gt. real(1E-07,TEST_AST_KIND_)) then
        print *, ref_fracday, " @ ", i, " does not equal ", fracday
        call mpp_error(FATAL, "test_diurnal_solar: 0d fracday value does not match reference value")
      end if

      if (abs(ref_rrsun - rrsun) .gt. real(1E-07,TEST_AST_KIND_)) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, "test_diurnal_solar: 0d rrsun value does not match reference value")
      end if 
    end do

  end subroutine test_diurnal_solar_0d

end program test_diurnal_solar_cal