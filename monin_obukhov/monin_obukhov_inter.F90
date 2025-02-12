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
!> @defgroup monin_obukhov_inter monin_obukhov_inter
!> @ingroup monin_obukhov
!> @brief Utility routines to be used in @ref monin_obukhov_mod

!> @addtogroup monin_obukhov_inter
!> @{
module monin_obukhov_inter

use platform_mod,    only: r4_kind, r8_kind
implicit none
private


public :: monin_obukhov_diff
public :: monin_obukhov_drag_1d
public :: monin_obukhov_solve_zeta
public :: monin_obukhov_derivative_t
public :: monin_obukhov_derivative_m
public :: monin_obukhov_profile_1d
public :: monin_obukhov_integral_m
public :: monin_obukhov_integral_tq
public :: monin_obukhov_stable_mix

interface monin_obukhov_diff
    module procedure monin_obukhov_diff_r4, monin_obukhov_diff_r8
end interface monin_obukhov_diff

interface monin_obukhov_drag_1d
    module procedure monin_obukhov_drag_1d_r4, monin_obukhov_drag_1d_r8
end interface monin_obukhov_drag_1d

interface monin_obukhov_solve_zeta
    module procedure monin_obukhov_solve_zeta_r4, monin_obukhov_solve_zeta_r8
end interface monin_obukhov_solve_zeta

interface monin_obukhov_derivative_t
    module procedure monin_obukhov_derivative_t_r4, monin_obukhov_derivative_t_r8
end interface monin_obukhov_derivative_t

interface monin_obukhov_derivative_m
    module procedure monin_obukhov_derivative_m_r4, monin_obukhov_derivative_m_r8
end interface monin_obukhov_derivative_m

interface monin_obukhov_profile_1d
    module procedure monin_obukhov_profile_1d_r4, monin_obukhov_profile_1d_r8
end interface monin_obukhov_profile_1d

interface monin_obukhov_integral_m
    module procedure monin_obukhov_integral_m_r4, monin_obukhov_integral_m_r8
end interface monin_obukhov_integral_m

interface monin_obukhov_integral_tq
    module procedure monin_obukhov_integral_tq_r4, monin_obukhov_integral_tq_r8
end interface monin_obukhov_integral_tq

interface monin_obukhov_stable_mix
    module procedure monin_obukhov_stable_mix_r4, monin_obukhov_stable_mix_r8
end interface monin_obukhov_stable_mix

contains

#include "monin_obukhov_inter_r4.fh"
#include "monin_obukhov_inter_r8.fh"

end module monin_obukhov_inter
!> @}
! close documentation grouping
