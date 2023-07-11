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
!> @brief Unit test for astronomy/annual_mean_solar_2level
!> @email gfdl.climate.model.info@noaa.gov
!> @description Performs calculations done in astronomy_mod using the daily_mean_solar
!> 2level subroutine using 32 and 64 bit reals

program test_annual_2level

  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only : mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use astronomy_mod,     only: astronomy_init, set_orbital_parameters, get_orbital_parameters, &
                               diurnal_solar, daily_mean_solar, annual_mean_solar
  use time_manager_mod,  only: JULIAN, set_calendar_type
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind
  use fms_string_utils_mod, only: stringify

  implicit none

  real(kind=r8_kind), parameter :: twopi = 2.0_r8_kind * real(PI, r8_kind)
  real(kind=r8_kind), parameter :: deg_to_rad  = twopi/360.0_r8_kind
  integer                       :: num_angles = 3600


  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call test_annual_mean_solar_2level

  call fms_end

  contains

  ! create reference values to test annual_mean_solar_2level
  subroutine create_ref_values_2level(jst, jnd, lat, cosz, solar, fracday, rrsun)

    implicit none
    integer, parameter                                   :: n = 16
    integer, intent(in)                                  :: jst, jnd
    real(kind=TEST_AST_KIND_), intent(in), dimension(n)  :: lat
    real(kind=TEST_AST_KIND_), intent(out), dimension(n) :: cosz, solar, fracday
    real(kind=TEST_AST_KIND_), intent(out)               :: rrsun
    real(kind=TEST_AST_KIND_)                            :: t
    real(kind=TEST_AST_KIND_), dimension(n)              :: s, z
    integer                                              :: i
    integer, parameter                                   :: lkind = TEST_AST_KIND_

    solar = 0.0_lkind
    cosz  = 0.0_lkind

    do i = 1, num_angles
      t = real((i-1),TEST_AST_KIND_) * real(twopi,TEST_AST_KIND_) / real(num_angles,TEST_AST_KIND_)
      call daily_mean_solar(lat, t, z, fracday, rrsun)
      s = z * rrsun * fracday
      solar = solar + s
      cosz = cosz + z * s
    end do

    solar = solar / num_angles
    cosz = cosz / num_angles

    where (solar .eq. 0.0_lkind)
      cosz    = 0.0_lkind
      fracday = 0.0_lkind
    elsewhere
      cosz    = cosz/solar
      fracday = solar/cosz
    end where

    rrsun = 1.0_lkind

  end subroutine create_ref_values_2level

  subroutine test_annual_mean_solar_2level

    implicit none
    integer                                  :: jst, jnd
    real(kind=TEST_AST_KIND_), dimension(16) :: lat
    real(kind=TEST_AST_KIND_), dimension(16) :: cosz, solar
    real(kind=TEST_AST_KIND_), dimension(16) :: ref_cosz, ref_solar, ref_fracday
    real(kind=TEST_AST_KIND_)                :: ref_rrsun
    integer                                  :: i
    integer, parameter                       :: lkind = TEST_AST_KIND_

    do i = 1, 16
      lat(i) = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    call create_ref_values_2level(jst, jnd, lat, ref_cosz, ref_solar, ref_fracday, ref_rrsun)
    call annual_mean_solar(lat, cosz, solar)

    do i = 1, 16

      if (abs(ref_cosz(i) - cosz(i)) .gt. real(1E-07,TEST_AST_KIND_)) then
        print *, ref_cosz(i), " @ ", i, " does not equal ", cosz(i)
        call mpp_error(FATAL, "test_annual_solar_2level: cosz value does not match referece value")
      end if

      if (abs(ref_solar(i) - solar(i)) .gt. real(1E-07,TEST_AST_KIND_)) then
        print *, ref_solar(i), " @ ", i, " does not equal ", solar(i)
        call mpp_error(FATAL, "test_annual_solar_2level: solar value does not match reference value")
      end if
    end do

  end subroutine test_annual_mean_solar_2level


end program test_annual_2level