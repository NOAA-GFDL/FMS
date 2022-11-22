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
!> @defgroup monin_obukhov_mod monin_obukhov_mod
!> @ingroup monin_obukhov
!> @brief Routines for computing surface drag coefficients
!! from data at the lowest model level
!! and for computing the profile of fields
!! between the lowest model level and the ground
!! using Monin-Obukhov scaling

module monin_obukhov_mod

use constants_mod, only: grav, vonkarm
use mpp_mod,       only: input_nml_file
use fms_mod,       only: error_mesg, FATAL, check_nml_error,   &
                         mpp_pe, mpp_root_pe, stdlog, &
                         write_version_number
use monin_obukhov_inter, only: monin_obukhov_diff, monin_obukhov_drag_1d, &
                               monin_obukhov_profile_1d, monin_obukhov_stable_mix
use platform_mod
implicit none
private

!=======================================================================
 public :: monin_obukhov_init
 public :: monin_obukhov_end
 public :: mo_drag
 public :: mo_profile
 public :: mo_diff
 public :: stable_mix
!=======================================================================

!> @brief Compute surface drag coefficients
!> @ingroup monin_obukhov_mod

interface monin_obukhov_init
    module procedure monin_obukhov_init_r4, monin_obukhov_init_r8
end interface monin_obukhov_init

interface mo_drag
    module procedure mo_drag_0d_r4, mo_drag_0d_r8
    module procedure mo_drag_1d_r4, mo_drag_1d_r8
    module procedure mo_drag_2d_r4, mo_drag_2d_r8
end interface mo_drag

!> @ingroup monin_obukhov_mod
interface mo_profile
    module procedure mo_profile_0d_r4, mo_profile_0d_r8
    module procedure mo_profile_1d_r4, mo_profile_1d_r8
    module procedure mo_profile_2d_r4, mo_profile_2d_r8 
    module procedure mo_profile_0d_n_r4, mo_profile_0d_n_r8
    module procedure mo_profile_1d_n_r4, mo_profile_1d_n_r8
    module procedure mo_profile_2d_n_r4, mo_profile_2d_n_r8
end interface mo_profile

!> @ingroup monin_obukhov_mod
interface mo_diff
    module procedure mo_diff_0d_n_r4, mo_diff_0d_n_r8
    module procedure mo_diff_0d_1_r4, mo_diff_0d_1_r8
    module procedure mo_diff_1d_n_r4, mo_diff_1d_n_r8
    module procedure mo_diff_1d_1_r4, mo_diff_1d_1_r8
    module procedure mo_diff_2d_n_r4, mo_diff_2d_n_r8
    module procedure mo_diff_2d_1_r4, mo_diff_2d_1_r8
end interface mo_diff

!> @ingroup monin_obukhov_mod
interface stable_mix
    module procedure stable_mix_0d_r4, stable_mix_0d_r8
    module procedure stable_mix_1d_r4, stable_mix_1d_r8
    module procedure stable_mix_2d_r4, stable_mix_2d_r8
    module procedure stable_mix_3d_r4, stable_mix_3d_r8
end interface stable_mix

!> @addtogroup monin_obukhov_mod
!> @{

!-----------------------------------------------------------------------
! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

!=======================================================================

!  DEFAULT VALUES OF NAMELIST PARAMETERS:

real(r8_kind)    :: rich_crit      = 2.0
real(r8_kind)    :: drag_min_heat  = 1.e-05
real(r8_kind)    :: drag_min_moist = 1.e-05
real(r8_kind)    :: drag_min_mom   = 1.e-05
logical          :: neutral        = .false.
integer          :: stable_option  = 1
real(r8_kind)    :: zeta_trans     = 0.5
logical          :: new_mo_option  = .false.


namelist /monin_obukhov_nml/ rich_crit, neutral, drag_min_heat, &
                             drag_min_moist, drag_min_mom,      &
                             stable_option, zeta_trans, new_mo_option !miz

!=======================================================================

!  MODULE VARIABLES

real(r8_kind), parameter  :: small  = 1.e-04
real(r8_kind)             :: b_stab, r_crit, lambda, rich_trans
real(r8_kind)             :: sqrt_drag_min_heat, sqrt_drag_min_moist, sqrt_drag_min_mom
logical                   :: module_is_initialized = .false.

contains

subroutine monin_obukhov_end

module_is_initialized = .false.

end subroutine monin_obukhov_end

!=======================================================================
#include<monin_obukhov.inc>

end module monin_obukhov_mod
!> @}
! close documentation grouping
