!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
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
