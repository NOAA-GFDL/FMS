! -*-f90-*-

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

!> @{

#undef MO_DRAG_1D_
#define MO_DRAG_1D_ mo_drag_1d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_PROFILE_1D_
#define MO_PROFILE_1D_ mo_profile_1d_r4
#undef FMS_MO_KIND_ 
#define FMS_MO_KIND_ r4_kind

#undef STABLE_MIX_3D_
#define STABLE_MIX_3D_ stable_mix_3d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DIFF_2D_N_
#define MO_DIFF_2D_N_ mo_diff_2d_n_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef SOLVE_ZETA_
#define SOLVE_ZETA_ solve_zeta_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DERIVATIVE_M_
#define MO_DERIVATIVE_M_ mo_derivative_m_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DERIVATIVE_T_
#define MO_DERIVATIVE_T_ mo_derivative_t_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_INTEGRAL_TQ_
#define MO_INTEGRAL_TQ_ mo_integral_tq_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_INTEGRAL_M_
#define MO_INTEGRAL_M_ mo_integral_m_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DRAG_2D_
#define MO_DRAG_2D_ mo_drag_2d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DRAG_0D_
#define MO_DRAG_0D_ mo_drag_0d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_PROFILE_2D_
#define MO_PROFILE_2D_ mo_profile_2d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_PROFILE_0D_
#define MO_PROFILE_0D_ mo_profile_0d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_PROFILE_1D_N_
#define MO_PROFILE_1D_N_ mo_profile_1d_n_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_PROFILE_0D_N_
#define MO_PROFILE_0D_N_ mo_profile_0d_n_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_PROFILE_2D_N_
#define MO_PROFILE_2D_N_ mo_profile_2d_n_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DIFF_2D_1_
#define MO_DIFF_2D_1_ mo_diff_2d_1_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DIFF_1D_1_
#define MO_DIFF_1D_1_ mo_diff_1d_1_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DIFF_1D_N_
#define MO_DIFF_1D_N_ mo_diff_1d_n_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DIFF_0D_1_
#define MO_DIFF_0D_1_ mo_diff_0d_1_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef MO_DIFF_0D_N_
#define MO_DIFF_0D_N_ mo_diff_0d_n_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef STABLE_MIX_2D_
#define STABLE_MIX_2D_ stable_mix_2d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef STABLE_MIX_1D_
#define STABLE_MIX_1D_ stable_mix_1d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind

#undef STABLE_MIX_0D_
#define STABLE_MIX_0D_ stable_mix_0d_r4
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r4_kind
#include <monin_obukhov.fh>

#undef MONIN_OBUKHOV_INIT_
#define MONIN_OBUKHOV_INIT_ monin_obukhov_init_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DRAG_1D_
#define MO_DRAG_1D_ mo_drag_1d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_PROFILE_1D_
#define MO_PROFILE_1D_ mo_profile_1d_r8
#undef FMS_MO_KIND_ 
#define FMS_MO_KIND_ r8_kind

#undef STABLE_MIX_3D_
#define STABLE_MIX_3D_ stable_mix_3d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DIFF_2D_N_
#define MO_DIFF_2D_N_ mo_diff_2d_n_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef SOLVE_ZETA_
#define SOLVE_ZETA_ solve_zeta_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DERIVATIVE_M_
#define MO_DERIVATIVE_M_ mo_derivative_m_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DERIVATIVE_T_
#define MO_DERIVATIVE_T_ mo_derivative_t_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_INTEGRAL_TQ_
#define MO_INTEGRAL_TQ_ mo_integral_tq_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_INTEGRAL_M_
#define MO_INTEGRAL_M_ mo_integral_m_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DRAG_2D_
#define MO_DRAG_2D_ mo_drag_2d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DRAG_0D_
#define MO_DRAG_0D_ mo_drag_0d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_PROFILE_2D_
#define MO_PROFILE_2D_ mo_profile_2d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_PROFILE_0D_
#define MO_PROFILE_0D_ mo_profile_0d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_PROFILE_1D_N_
#define MO_PROFILE_1D_N_ mo_profile_1d_n_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_PROFILE_0D_N_
#define MO_PROFILE_0D_N_ mo_profile_0d_n_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_PROFILE_2D_N_
#define MO_PROFILE_2D_N_ mo_profile_2d_n_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DIFF_2D_1_
#define MO_DIFF_2D_1_ mo_diff_2d_1_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DIFF_1D_1_
#define MO_DIFF_1D_1_ mo_diff_1d_1_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DIFF_1D_N_
#define MO_DIFF_1D_N_ mo_diff_1d_n_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DIFF_0D_1_
#define MO_DIFF_0D_1_ mo_diff_0d_1_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef MO_DIFF_0D_N_
#define MO_DIFF_0D_N_ mo_diff_0d_n_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef STABLE_MIX_2D_
#define STABLE_MIX_2D_ stable_mix_2d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef STABLE_MIX_1D_
#define STABLE_MIX_1D_ stable_mix_1d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef STABLE_MIX_0D_
#define STABLE_MIX_0D_ stable_mix_0d_r8
#undef FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind
#include <monin_obukhov.fh>

!> @}
! close documentation grouping