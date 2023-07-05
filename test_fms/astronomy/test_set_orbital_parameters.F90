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
!> @description Tests that astronomy_mod breaks when orbital 
!> parameters are out of bounds using 32 and 64 bit reals

program test_set_parameters

    use fms_mod,           only: fms_init, fms_end
    use mpp_mod,           only : mpp_error, FATAL, stdout, mpp_init, mpp_exit
    use astronomy_mod,     only: astronomy_init, set_orbital_parameters, get_orbital_parameters, &
                                 diurnal_solar, daily_mean_solar, annual_mean_solar
    use time_manager_mod,  only: JULIAN, set_calendar_type
    use platform_mod,      only: r4_kind, r8_kind

    implicit none

    real(kind=TEST_AST_KIND_) :: ecc   = 0.01671
    real(kind=TEST_AST_KIND_) :: obliq = 23.439
    real(kind=TEST_AST_KIND_) :: per   = 102.932


    call fms_init()
    call set_calendar_type(JULIAN)
    call astronomy_init

    call test_set_orbital_parameters

    call fms_end

    contains

    !Test set_orbital_parameters breaks when arguments are not in bounds (XFAIL)
    subroutine test_set_orbital_parameters

        implicit none
        real(kind=TEST_AST_KIND_) :: ecc_in_check
        real(kind=TEST_AST_KIND_) :: obliq_in_check
        real(kind=TEST_AST_KIND_) :: per_in_check
        integer, parameter :: lkind = TEST_AST_KIND_

        ! Eccentricity of Earth's orbit not in between 0.0 and 0.99
        ecc_in_check = -1.0_lkind
        call set_orbital_parameters(ecc_in_check, obliq, per)
        ecc_in_check = 1.0_lkind
        call set_orbital_parameters(ecc_in_check, obliq, per)

        ! Obliquity not in between -90 and 90
        obliq_in_check = -91.0_lkind
        call set_orbital_parameters(ecc, obliq_in_check, per)
        obliq_in_check = 91.0_lkind
        call set_orbital_parameters(ecc, obliq_in_check, per)

        ! Perihelion not in between 0.0 and 360.0
        per_in_check = -1.0_lkind
        call set_orbital_parameters(ecc, obliq, per_in_check)
        per_in_check = 361.0_lkind
        call set_orbital_parameters(ecc, obliq, per_in_check)


    end subroutine test_set_orbital_parameters

end program test_set_parameters