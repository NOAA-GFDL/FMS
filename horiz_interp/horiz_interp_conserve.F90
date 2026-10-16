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
!> @defgroup horiz_interp_conserve_mod horiz_interp_conserve_mod
!! @ingroup horiz_interp
!! @{
!!
!! @author Bruce Wyman, Zhi Liang
!!
!! @parblock
!! Horiz_interp_conserve_mod contains methods called from horiz_interp_mod to
!! interpolate data on any logically rectangular grid to any rectangular grid
!! with order 1 conservative method.  Users are recommended to use the top-level
!! module horiz_interp_mod with "interp" set to "conservative for conservative
!! interpolation.
!! @endparblock

module horiz_interp_conserve_mod

  use platform_mod,          only: r4_kind, r8_kind
  use mpp_mod,               only: mpp_send, mpp_recv, mpp_pe, mpp_root_pe, mpp_npes
  use mpp_mod,               only: mpp_error, FATAL,  mpp_sync_self
  use mpp_mod,               only: COMM_TAG_1, COMM_TAG_2
  use fms_mod,               only: write_version_number
  use grid2_mod,             only: get_great_circle_algorithm
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type, CONSERVE


  implicit none
  private


  !> Generic interface to compute interpolation weights and mapping indices.
  !! Member subroutines and their main arguments are:
  !! horiz_interp_conserve_new_1dx1d_r4:
  !!  input and output grids are provided as 1D arrays in 32-bit precision.
  !! horiz_interp_conserve_new_1dx1d_r8:
  !!  input and output grids are provided as 1D arrays in 64-bit precision.
  !! horiz_interp_conserve_new_2dx1d_r4:
  !!  input grid is provided as 2D arrays, output grid as 1D arrays, both in 32-bit precision.
  !! horiz_interp_conserve_new_2dx1d_r8:
  !!  input grid is provided as 2D arrays, output grid as 1D arrays, both in 64-bit precision.
  !! horiz_interp_conserve_new_1dx2d_r4:
  !!   input grid is provided as 1D arrays, output grid as 2D arrays, both in 32-bit precision.
  !! horiz_interp_conserve_new_1dx2d_r8:
  !!   input grid is provided as 1D arrays, output grid as 2D arrays, both in 64-bit precision.
  !! horiz_interp_conserve_new_2dx2d_r4:
  !!   input and output grids are provided as 2D arrays in 32-bit precision.
  !! horiz_interp_conserve_new_2dx2d_r8:
  !!   input and output grids are provided as 2D arrays in 64-bit precision.
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

  !> Generic interface to interpolate data from source grid to target grid from the weights
  !! stored in Interp.  Calls horiz_interp_conserve_r4 if input and output data are 2D arrays
  !! in 32-bit precision.  Calls horiz_interp_conserve_r8 if input and output data are 2D
  !! arrays in 64-bit precision.  Horiz_interp_conserve_r* will subsequently call
  !! horiz_interp_conserve_version1 or horiz_interp_conserve_version2 depending on value of
  !! Interp%version.
  interface horiz_interp_conserve
    module procedure horiz_interp_conserve_r4
    module procedure horiz_interp_conserve_r8
  end interface

  !> Internally used generic interface when interpolating data with Interp
  !! generated from horiz_interp_new_1dx1d_r*.  Not frequently used.
  interface data_sum
    module procedure data_sum_r4
    module procedure data_sum_r8
  end interface

  !> Internally used generic interface when interpolating data with Interp
  !! generated from horiz_interp_new_1dx1d_r*.  Not frequently used.
  interface stats
    module procedure stats_r4
    module procedure stats_r8
  end interface

  !> Generic inteface called from horiz_interp_conserve to interpolate
  !! data using Interp populated from horiz_interp_new_1dx1d_r*.
  interface horiz_interp_conserve_version1
    module procedure horiz_interp_conserve_version1_r8
    module procedure horiz_interp_conserve_version1_r4
  end interface

  !> Generic interface called from horiz_interp_conserve to interpolate
  !! 2d data.  Calls horiz_interp_conserve_version2_r4 when data are 2D
  !! arrays in 32-bit precision.  Calls horiz_interp_conserve_version2_r8 when
  !! data are 2D arrays in 64-bit precision.
  interface horiz_interp_conserve_version2
    module procedure horiz_interp_conserve_version2_r8
    module procedure horiz_interp_conserve_version2_r4
  end interface

  public :: horiz_interp_conserve_init
  public :: horiz_interp_conserve_new, horiz_interp_conserve, horiz_interp_conserve_del

  integer :: pe !< is the current pe
  integer :: root_pe !< is the root_pe

#include<file_version.h>

  !! TODO:  how to set great_circle_algorithm
  logical :: module_is_initialized = .FALSE. !< is a flag to prevent module re-initialization.
  logical :: great_circle_algorithm = .false. !< turns on the great_circle_algorithm to compute grid cell areas.

contains

  !> @parblock
  !! Initializes horiz_interp_conserve_mod.  Called from horiz_interp_init in horiz_interp_mod.
  !! @endparblock
  subroutine horiz_interp_conserve_init

    if(module_is_initialized) return
    call write_version_number("HORIZ_INTERP_CONSERVE_MOD", version)

    great_circle_algorithm = get_great_circle_algorithm()

    module_is_initialized = .true.

  end subroutine horiz_interp_conserve_init

  !> @parblock
  !! Deallocates arrays holding conservative order 1 interpolation weights and mapping indices
  !! in Interp.  Resets %is_allocated to .false.  Called from horiz_interp_del in horiz_interp_mod.
  !! @endparblock
  subroutine horiz_interp_conserve_del ( Interp )

    type (horiz_interp_type), intent(inout) :: Interp !< will be reset with deallocated memory

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
