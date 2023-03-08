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

#undef FMS_DIAG_FBU_DATA_TYPE_
#define FMS_DIAG_FBU_DATA_TYPE_ REAL(r4_kind)
#undef FMS_DIAG_FBU_PNAME_
#define FMS_DIAG_FBU_PNAME_ fieldbuff_update_r4
#undef  FMS_DIAG_FBU_3D_PNAME_
#define FMS_DIAG_FBU_3D_PNAME_ fieldbuff_update_3d_r4
#undef  FMS_DIAG_FBCF_PNAME_
#define FMS_DIAG_FBCF_PNAME_ fieldbuff_copy_fieldvals_r4
#undef  FMS_DIAG_FBCF_3D_PNAME_
#define FMS_DIAG_FBCF_3D_PNAME_ fieldbuff_copy_fieldvals_3d_r4
#undef  FMS_DIAG_FBCM_PNAME_
#define FMS_DIAG_FBCM_PNAME_ fieldbuff_copy_missvals_r4
#undef  FMS_DIAG_FBCM_3D_PNAME_
#define FMS_DIAG_FBCM_3D_PNAME_ fieldbuff_copy_missvals_3d_r4
#include <fms_diag_fieldbuff_update.fh>

#undef FMS_DIAG_FBU_DATA_TYPE_
#define FMS_DIAG_FBU_DATA_TYPE_ REAL(r8_kind)
#undef FMS_DIAG_FBU_PNAME_
#define FMS_DIAG_FBU_PNAME_ fieldbuff_update_r8
#undef  FMS_DIAG_FBU_3D_PNAME_
#define FMS_DIAG_FBU_3D_PNAME_ fieldbuff_update_3d_r8
#undef  FMS_DIAG_FBCF_PNAME_
#define FMS_DIAG_FBCF_PNAME_ fieldbuff_copy_fieldvals_r8
#undef  FMS_DIAG_FBCF_3D_PNAME_
#define FMS_DIAG_FBCF_3D_PNAME_ fieldbuff_copy_fieldvals_3d_r8
#undef  FMS_DIAG_FBCM_PNAME_
#define FMS_DIAG_FBCM_PNAME_ fieldbuff_copy_missvals_r8
#undef  FMS_DIAG_FBCM_3D_PNAME_
#define FMS_DIAG_FBCM_3D_PNAME_ fieldbuff_copy_missvals_3d_r8
#include <fms_diag_fieldbuff_update.fh>
