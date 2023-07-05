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
!> 2d subroutine using 32 and 64 bit reals

program test_annual_2d

  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only : mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use astronomy_mod,     only: astronomy_init, set_orbital_parameters, get_orbital_parameters, &
                               diurnal_solar, daily_mean_solar, annual_mean_solar
  use time_manager_mod,  only: JULIAN, set_calendar_type
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind
  use fms_string_utils_mod, only: stringify

  implicit none

  real(kind=TEST_AST_KIND_), parameter :: twopi = real(2.0,TEST_AST_KIND_) * real(PI, TEST_AST_KIND_)
  real(kind=TEST_AST_KIND_), parameter :: deg_to_rad  = twopi/real(360.0,TEST_AST_KIND_)
  integer                  :: num_angles = 3600


  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call test_annual_mean_solar_2d

  call fms_end

  contains

  ! create reference values to test annual_mean_solar_2d
  subroutine create_ref_values_2d(jst, jnd, lat, cosz, solar, fracday, rrsun)

    implicit none
    integer, intent(in)                                    :: jst, jnd
    real(kind=TEST_AST_KIND_), intent(in), dimension(4,4)  :: lat
    real(kind=TEST_AST_KIND_), intent(out), dimension(4,4) :: cosz, solar, fracday
    real(kind=TEST_AST_KIND_), intent(out)                 :: rrsun
    real(kind=TEST_AST_KIND_)                              :: t
    real(kind=TEST_AST_KIND_), dimension(4,4)              :: s, z
    integer                                                :: n
    integer, parameter                                     :: lkind = TEST_AST_KIND_

    solar = 0.0_lkind
    cosz = 0.0_lkind

    do n = 1, num_angles
      t = real((n-1),TEST_AST_KIND_) * real(twopi,TEST_AST_KIND_) / real(num_angles,TEST_AST_KIND_)
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

  end subroutine create_ref_values_2d
  
  subroutine test_annual_mean_solar_2d

    implicit none
    integer                                   :: jst, jnd
    real(kind=TEST_AST_KIND_), dimension(4,4) :: lat
    real(kind=TEST_AST_KIND_), dimension(4,4) :: cosz, solar, fracday
    real(kind=TEST_AST_KIND_), dimension(4,4) :: ref_cosz, ref_solar, ref_fracday
    real(kind=TEST_AST_KIND_)                 :: ref_rrsun, rrsun
    integer                                   :: i, j, counter
    integer, parameter                        :: lkind = TEST_AST_KIND_

    jst = 0.0_lkind
    jnd = 360.0_lkind * deg_to_rad

    do i = 1, 4
      do j = 1, 4
        lat(j,i) = real(counter, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
        counter = counter + 1
      end do
    end do

    call create_ref_values_2d(jst, jnd, lat, ref_cosz, ref_solar, ref_fracday, ref_rrsun)
    call annual_mean_solar(jst, jnd, lat, cosz, solar, fracday, rrsun)

    do i = 1, 4
      do j = 1, 4

        if (abs(ref_cosz(i,j) - cosz(i,j)) .gt. real(1E-07,TEST_AST_KIND_)) then
          print *, ref_cosz(i,j), " @ ", i,j, " does not equal ", cosz(i,j)
          call mpp_error(FATAL, "test_annual_solar_2d: cosz value does not match reference value")
        end if

        if (abs(ref_solar(i,j) - solar(i,j)) .gt. real(1E-07,TEST_AST_KIND_)) then
          print *, ref_solar(i,j), " @ ", i,j, " does not equal ", solar(i,j)
          call mpp_error(FATAL, "test_annual_solar_2d: solar value does not match reference value")
        end if

        if (abs(ref_fracday(i,j) - fracday(i,j)) .gt. real(1E-07,TEST_AST_KIND_)) then
          print *, ref_fracday(i,j), " @ ", i,j, " does not equal ", fracday(i,j)
          call mpp_error(FATAL, "test_annual_solar_2d: fracday value does not match reference value")
        end if
      end do
    end do

    if (abs(ref_rrsun - rrsun) .gt. real(1E-07,TEST_AST_KIND_)) then
      print *, ref_solar, " does not equal ", rrsun
      call mpp_error(FATAL, "test_annual_solar_2d: rrsun value does not match reference value")
    end if
    
  end subroutine test_annual_mean_solar_2d


end program test_annual_2d