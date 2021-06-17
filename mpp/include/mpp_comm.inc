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
!> @file
!> @ingroup mpp
!> @brief Imports checksum, gather, and scatter routines used for communications and calculations
!! between PE's in @ref mpp_mod

#if defined(use_libMPI)
#include <mpp_comm_mpi.inc>
#else
#include <mpp_comm_nocomm.inc>
#endif

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_1d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i8_1d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i8_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_2d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i8_2d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i8_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_3d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i8_3d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i8_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_4d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i8_4d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i8_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_5d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i8_5d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i8_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_1d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i4_1d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_2d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i4_2d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_3d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i4_3d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_4d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i4_4d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_5d
#undef MPP_CHKSUM_INT_RMASK_
#define MPP_CHKSUM_INT_RMASK_ mpp_chksum_i4_5d_rmask
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_0d
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_1d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_2d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_3d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_4d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_5d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#ifdef OVERLOAD_C8
#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_0d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_1d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_2d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_3d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_4d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_5d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i8_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>
#endif

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_0d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_1d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_2d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_3d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_4d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_5d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#ifdef OVERLOAD_C4
#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_0d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_1d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_2d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_3d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_4d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_5d
#undef MPP_TRANSFER_KIND_
#define MPP_TRANSFER_KIND_ i4_kind
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>
#endif

!#################################################
#undef MPP_GATHER_1D_
#undef MPP_GATHER_1DV_
#undef MPP_TYPE_
#define MPP_TYPE_ logical
#define MPP_GATHER_1D_ mpp_gather_logical_1d
#define MPP_GATHER_1DV_ mpp_gather_logical_1dv
#undef  MPP_GATHER_PELIST_2D_
#undef  MPP_GATHER_PELIST_3D_
#define MPP_GATHER_PELIST_2D_ mpp_gather_pelist_logical_2d
#define MPP_GATHER_PELIST_3D_ mpp_gather_pelist_logical_3d
#include <mpp_gather.h>

#undef MPP_GATHER_1D_
#undef MPP_GATHER_1DV_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#define MPP_GATHER_1D_ mpp_gather_int4_1d
#define MPP_GATHER_1DV_ mpp_gather_int4_1dv
#undef  MPP_GATHER_PELIST_2D_
#undef  MPP_GATHER_PELIST_3D_
#define MPP_GATHER_PELIST_2D_ mpp_gather_pelist_int4_2d
#define MPP_GATHER_PELIST_3D_ mpp_gather_pelist_int4_3d
#include <mpp_gather.h>

#undef MPP_GATHER_1D_
#undef MPP_GATHER_1DV_
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#define MPP_GATHER_1D_ mpp_gather_real4_1d
#define MPP_GATHER_1DV_ mpp_gather_real4_1dv
#undef  MPP_GATHER_PELIST_2D_
#undef  MPP_GATHER_PELIST_3D_
#define MPP_GATHER_PELIST_2D_ mpp_gather_pelist_real4_2d
#define MPP_GATHER_PELIST_3D_ mpp_gather_pelist_real4_3d
#include <mpp_gather.h>

#undef MPP_GATHER_1D_
#undef MPP_GATHER_1DV_
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#define MPP_GATHER_1D_ mpp_gather_real8_1d
#define MPP_GATHER_1DV_ mpp_gather_real8_1dv
#undef  MPP_GATHER_PELIST_2D_
#undef  MPP_GATHER_PELIST_3D_
#define MPP_GATHER_PELIST_2D_ mpp_gather_pelist_real8_2d
#define MPP_GATHER_PELIST_3D_ mpp_gather_pelist_real8_3d
#include <mpp_gather.h>

!#################################################
#undef  MPP_SCATTER_PELIST_2D_
#undef  MPP_SCATTER_PELIST_3D_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#define MPP_SCATTER_PELIST_2D_ mpp_scatter_pelist_int4_2d
#define MPP_SCATTER_PELIST_3D_ mpp_scatter_pelist_int4_3d
#include <mpp_scatter.h>

#undef  MPP_SCATTER_PELIST_2D_
#undef  MPP_SCATTER_PELIST_3D_
#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#define MPP_SCATTER_PELIST_2D_ mpp_scatter_pelist_real4_2d
#define MPP_SCATTER_PELIST_3D_ mpp_scatter_pelist_real4_3d
#include <mpp_scatter.h>

#undef  MPP_SCATTER_PELIST_2D_
#undef  MPP_SCATTER_PELIST_3D_
#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#define MPP_SCATTER_PELIST_2D_ mpp_scatter_pelist_real8_2d
#define MPP_SCATTER_PELIST_3D_ mpp_scatter_pelist_real8_3d
#include <mpp_scatter.h>

