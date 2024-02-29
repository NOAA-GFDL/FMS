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
!> @defgroup horiz_interp_conserve_mod horiz_interp_conserve_mod
!> @ingroup horiz_interp
!> @brief Performs spatial interpolation between grids using conservative interpolation
!!
!> @author Bruce Wyman, Zhi Liang
!!
!> This module can conservatively interpolate data from any logically rectangular grid
!! to any rectangular grid. The interpolation scheme is area-averaging
!! conservative scheme. There is an optional mask field for missing input data in both
!! horiz_interp__conserveinit and horiz_interp_conserve. For efficiency purpose, mask should only be
!! kept in horiz_interp_init (will remove the mask in horiz_interp in the future).
!! There are 1-D and 2-D version of horiz_interp_conserve_init for 1-D and 2-D grid.
!! There is a optional argument mask in horiz_interp_conserve_init_2d and no mask should
!! to passed into horiz_interp_conserv. optional argument mask will not be passed into
!! horiz_interp_conserve_init_1d and optional argument mask may be passed into
!! horiz_interp_conserve (For the purpose of reproduce Memphis??? results).
!! An optional output mask field may be used in conjunction with the input mask to show
!! where output data exists.

module horiz_interp_conserve_mod

  use platform_mod,          only: r4_kind, r8_kind
  use mpp_mod,               only: mpp_send, mpp_recv, mpp_pe, mpp_root_pe, mpp_npes
  use mpp_mod,               only: mpp_error, FATAL,  mpp_sync_self
  use mpp_mod,               only: COMM_TAG_1, COMM_TAG_2
  use fms_mod,               only: write_version_number
  use grid2_mod,             only: get_great_circle_algorithm
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type


  implicit none
  private

  ! public interface


  !> @brief Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights.
  !!
  !> Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights
  !! for improved performance of multiple interpolations between
  !! the same grids.
  !! @param lon_in
  !!      Longitude (in radians) for source data grid.
  !!
  !! @param lat_in
  !!      Latitude (in radians) for source data grid.
  !!
  !! @param lon_out
  !!      Longitude (in radians) for destination data grid.
  !!
  !! @param lat_out
  !!      Latitude (in radians) for destination data grid.
  !!
  !! @param verbose
  !!      flag for the amount of print output.
  !!
  !! @param mask_in
  !!      Input mask.  must be the size (size(lon_in)-1, size(lon. The real value of
  !!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points
  !!      that should not be used or have missing data.
  !!
  !! @param mask_out
  !!      Output mask that specifies whether data was computed.
  !!
  !! @param Interp
  !!      A derived-type variable containing indices and weights used for subsequent
  !!      interpolations. To reinitialize this variable for a different grid-to-grid
  !!      interpolation you must first use the "horiz_interp_del" interface.
  !!
  !> @ingroup horiz_interp_conserve_mod
  interface horiz_interp_conserve_new
     module procedure horiz_interp_conserve_new_1dx1d_r4
     module procedure horiz_interp_conserve_new_1dx2d_r4
     module procedure horiz_interp_conserve_new_2dx1d_r4
     module procedure horiz_interp_conserve_new_2dx2d_r4
     module procedure horiz_interp_conserve_new_1dx1d_r8
     module procedure horiz_interp_conserve_new_1dx2d_r8
     module procedure horiz_interp_conserve_new_2dx1d_r8
     module procedure horiz_interp_conserve_new_2dx2d_r8
  end interface

  interface horiz_interp_conserve
    module procedure horiz_interp_conserve_r4
    module procedure horiz_interp_conserve_r8
  end interface

!> private helper routines
  interface data_sum
    module procedure data_sum_r4
    module procedure data_sum_r8
  end interface

  interface stats
    module procedure stats_r4
    module procedure stats_r8
  end interface

  interface horiz_interp_conserve_version1
    module procedure horiz_interp_conserve_version1_r8
    module procedure horiz_interp_conserve_version1_r4
  end interface

  interface horiz_interp_conserve_version2
    module procedure horiz_interp_conserve_version2_r8
    module procedure horiz_interp_conserve_version2_r4
  end interface



  !> @addtogroup horiz_interp_conserve_mod
  !> @{
  public :: horiz_interp_conserve_init
  public :: horiz_interp_conserve_new, horiz_interp_conserve, horiz_interp_conserve_del

  integer :: pe, root_pe
  !-----------------------------------------------------------------------
  ! Include variable "version" to be written to log file.
#include<file_version.h>
  logical            :: module_is_initialized = .FALSE.

  logical         :: great_circle_algorithm = .false.

contains

  !> Writes version number to logfile.
  subroutine horiz_interp_conserve_init

    if(module_is_initialized) return
    call write_version_number("HORIZ_INTERP_CONSERVE_MOD", version)

    great_circle_algorithm = get_great_circle_algorithm()

    module_is_initialized = .true.

  end subroutine horiz_interp_conserve_init

  !> Deallocates memory used by "HI_KIND_TYPE" variables.
  !! Must be called before reinitializing with horiz_interp_new.
  subroutine horiz_interp_conserve_del ( Interp )

    type (horiz_interp_type), intent(inout) :: Interp !< A derived-type variable returned by
                         !! previous call to horiz_interp_new. The input variable must have
                         !! allocated arrays. The returned variable will contain deallocated arrays.

    select case(Interp%version)
    case (1)
      if( Interp%horizInterpReals8_type%is_allocated) then
        if(allocated(Interp%horizInterpReals8_type%area_src)) deallocate(Interp%horizInterpReals8_type%area_src)
        if(allocated(Interp%horizInterpReals8_type%area_dst)) deallocate(Interp%horizInterpReals8_type%area_dst)
        if(allocated(Interp%horizInterpReals8_type%facj))     deallocate(Interp%horizInterpReals8_type%facj)
        if(allocated(Interp%jlat))                 deallocate(Interp%jlat)
        if(allocated(Interp%horizInterpReals8_type%faci))     deallocate(Interp%horizInterpReals8_type%faci)
        if(allocated(Interp%ilon))                 deallocate(Interp%ilon)
      else if( Interp%horizInterpReals4_type%is_allocated) then
        if(allocated(Interp%horizInterpReals4_type%area_src)) deallocate(Interp%horizInterpReals4_type%area_src)
        if(allocated(Interp%horizInterpReals4_type%area_dst)) deallocate(Interp%horizInterpReals4_type%area_dst)
        if(allocated(Interp%horizInterpReals4_type%facj))     deallocate(Interp%horizInterpReals4_type%facj)
        if(allocated(Interp%jlat))                 deallocate(Interp%jlat)
        if(allocated(Interp%horizInterpReals4_type%faci))     deallocate(Interp%horizInterpReals4_type%faci)
        if(allocated(Interp%ilon))                 deallocate(Interp%ilon)
      endif
    case (2)
      if( Interp%horizInterpReals8_type%is_allocated) then
        if(allocated(Interp%i_src)) deallocate(Interp%i_src)
        if(allocated(Interp%j_src)) deallocate(Interp%j_src)
        if(allocated(Interp%i_dst)) deallocate(Interp%i_dst)
        if(allocated(Interp%j_dst)) deallocate(Interp%j_dst)
        if(allocated(Interp%horizInterpReals8_type%area_frac_dst)) &
            deallocate(Interp%horizInterpReals8_type%area_frac_dst)
      else if( Interp%horizInterpReals4_type%is_allocated ) then
        if(allocated(Interp%i_src)) deallocate(Interp%i_src)
        if(allocated(Interp%j_src)) deallocate(Interp%j_src)
        if(allocated(Interp%i_dst)) deallocate(Interp%i_dst)
        if(allocated(Interp%j_dst)) deallocate(Interp%j_dst)
        if(allocated(Interp%horizInterpReals4_type%area_frac_dst)) &
            deallocate(Interp%horizInterpReals4_type%area_frac_dst)
       endif
    end select
    Interp%horizInterpReals4_type%is_allocated = .false.
    Interp%horizInterpReals8_type%is_allocated = .false.

  end subroutine horiz_interp_conserve_del

#include "horiz_interp_conserve_r4.fh"
#include "horiz_interp_conserve_r8.fh"

end module horiz_interp_conserve_mod
!> @}
! close documentation grouping
