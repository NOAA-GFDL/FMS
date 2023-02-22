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

!> @defgroup fms_diag_fieldbuff_update_mod fms_diag_fieldbuff_update_mod
!> @ingroup diag_manager
!> @brief fms_diag_fieldbuff_update_mod Contains routines for updating the
!! buffer (array) of field data statistics (e.g. average, rms) with new field data.
!!
!> @author Miguel Zuniga
!!
!! <TT>fms_diag_fieldbuff_update_mod</TT> contains routines for updating the buffer
!!(array) of field data statistics (e.g. average, rms) with new field data. These
!! routines are called by the send_data routines in the diag_manager.
!!
!> @file
!> @brief File for @ref fms_diag_fieldbuff_update_mod
!> @addtogroup fms_diag_fieldbuff_update_mod
!> @{
MODULE fms_diag_fieldbuff_update_mod
   USE platform_mod
   USE mpp_mod, ONLY: mpp_pe, mpp_root_pe
   USE time_manager_mod, ONLY: time_type
   USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdout, stdlog, write_version_number,fms_error_handler
   USE diag_data_mod, ONLY:  debug_diag_manager
   USE fms_diag_outfield_mod, ONLY: fmsDiagOutfieldIndex_type, fmsDiagOutfield_type
   USE diag_util_mod, ONLY: fms_diag_check_out_of_bounds
   USE fms_diag_time_reduction_mod, ONLY: fmsDiagTimeReduction_type
   USE fms_diag_elem_weight_procs_mod, ONLY: addwf
   USE fms_diag_bbox_mod, ONLY: fmsDiagIbounds_type

   implicit none

   !!TODO: (MDM) Remove commented integer versions.

   !> @brief Interface fieldbuff_update updates elements of field output buffer based on input field
   !! data and mathematical operations on the field data.
   !> @ingroup fms_diag_fieldbuff_update_mod
   interface fieldbuff_update
      !< r4 version of the interface
      module procedure fieldbuff_update_r4
      !< r8 version of the interface
      module procedure fieldbuff_update_r8
      !< r4 version of the interface, where the field is 3D
      module procedure fieldbuff_update_3d_r4
      !< r8 version of the interface
      module procedure fieldbuff_update_3d_r8
      !< i4 version of the interface, , where the field is 3D
      !module procedure fieldbuff_update_i4
      !< i8 version of the interface
     ! module procedure fieldbuff_update_i8
   end interface

   !> @brief Interface fieldbuff_copy_missvals updates elements of the field output buffer with
   !! the missvalue input argument.
   !> @ingroup fms_diag_fieldbuff_update_mod
   interface fieldbuff_copy_missvals
      !< r4 version of the interface
      module procedure fieldbuff_copy_missvals_r4
      !< r8 version of the interface
      module procedure fieldbuff_copy_missvals_r8
      !< r4 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_missvals_3d_r4
      !< r8 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_missvals_3d_r8
      !< i4 version of the interface
      !module procedure fieldbuff_copy_missvals_i4
      !< i8 version of the interface
      !module procedure fieldbuff_copy_missvals_i8
   end interface

   !> @brief Interface fieldbuff_copy_fieldvals updates elements of the field output buffer with
   !! copies of corresponding element values in the input field data.
   !> @ingroup fms_diag_fieldbuff_update_mod
   interface fieldbuff_copy_fieldvals
      !< r4 version of the interface
      module procedure fieldbuff_copy_fieldvals_r4
      !< r8 version of the interface
      module procedure fieldbuff_copy_fieldvals_r8
      !< r4 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_fieldvals_3d_r4
      !< r8 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_fieldvals_3d_r8
      !< i4 version of the interface
      !module procedure fieldbuff_copy_fieldvals_i4
      !< i8 version of the interface
      !module procedure fieldbuff_copy_fieldvals_i8
  end interface
contains

#include "fms_diag_fieldbuff_update.inc"

END MODULE fms_diag_fieldbuff_update_mod
!> @}
! close documentation grouping
