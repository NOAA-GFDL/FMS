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
!> @brief Unit test for astronomy/test_orbital_parameters
!> @email gfdl.climate.model.info@noaa.gov
!> @description Tests that astronomical parameters retrieves the orbital
!> parameters correctly using 32 and 64 bit reals


program test_get_parameters

  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only : mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use astronomy_mod,     only: astronomy_init, set_orbital_parameters, get_orbital_parameters, &
                               diurnal_solar, daily_mean_solar, annual_mean_solar
  use time_manager_mod,  only: JULIAN, set_calendar_type
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  real(kind=r8_kind) :: ecc   = 0.01671_r8_kind
  real(kind=r8_kind) :: obliq = 23.439_r8_kind
  real(kind=r8_kind) :: per   = 102.932_r8_kind


  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call test_get_orbital_parameters

  call fms_end

  contains

  subroutine test_get_orbital_parameters

    implicit none
    real(kind=r8_kind) :: ecc_check, obliq_check, per_check

    call get_orbital_parameters(ecc_check, obliq_check, per_check)

    if (ecc_check .ne. ecc) then
      call mpp_error(FATAL, "test_get_orbital_parameters: eccentricity precision was lost")
    else if (obliq_check .ne. obliq) then
      call mpp_error(FATAL, "test_get_orbital_parameters: obliquity precision was lost")
    else if (per .ne. per_check) then
      call mpp_error(FATAL, "test_get_orbital_parameters: perihelion precision was lost")
    else
    end if

  end subroutine test_get_orbital_parameters

end program test_get_parameters