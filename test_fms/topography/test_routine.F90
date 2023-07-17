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
!> @brief Unit test for astronomy/test_topography
!> @email gfdl.climate.model.info@noaa.gov
!> @description FILL THIS OUT

program test_top

  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only : mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use astronomy_mod,     only: astronomy_init, set_orbital_parameters, get_orbital_parameters, &
                               diurnal_solar, daily_mean_solar, annual_mean_solar
  use time_manager_mod,  only: JULIAN, set_calendar_type
  use platform_mod,      only: r4_kind, r8_kind

  implicit none


  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call test

  call fms_end

  contains

  subroutine test



  end subroutine test

end program test_top