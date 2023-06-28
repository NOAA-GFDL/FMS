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
!> @defgroup horiz_interp_type_mod horiz_interp_type_mod
!> @ingroup horiz_interp
!> @brief define derived data type that contains indices and weights used for subsequent
!! interpolations.
!> @author Zhi Liang

!> @addtogroup
!> @{
module horiz_interp_type_mod

use mpp_mod, only : mpp_send, mpp_recv, mpp_sync_self, mpp_error, FATAL
use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes
use mpp_mod, only : COMM_TAG_1, COMM_TAG_2
use platform_mod, only: r4_kind, r8_kind

implicit none
private


! parameter to determine interpolation method
 integer, parameter :: CONSERVE = 1
 integer, parameter :: BILINEAR = 2
 integer, parameter :: SPHERICA = 3
 integer, parameter :: BICUBIC  = 4

public :: CONSERVE, BILINEAR, SPHERICA, BICUBIC
public :: horiz_interp_type, stats, assignment(=)

!> @}

!> @ingroup horiz_interp_type_mod
interface assignment(=)
  module procedure horiz_interp_type_eq
end interface

!> @ingroup horiz_interp_type_mod
interface stats
  module procedure stats_r4
  module procedure stats_r8
end interface


!> real(8) pointers for use in horiz_interp_type
type horizInterpReals8_type
   real(kind=r8_kind),    dimension(:,:), allocatable   :: faci     !< weights for conservative scheme
   real(kind=r8_kind),    dimension(:,:), allocatable   :: facj     !< weights for conservative scheme
   real(kind=r8_kind),    dimension(:,:), allocatable   :: area_src !< area of the source grid
   real(kind=r8_kind),    dimension(:,:), allocatable   :: area_dst !< area of the destination grid
   real(kind=r8_kind),    dimension(:,:,:), allocatable :: wti      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r8_kind),    dimension(:,:,:), allocatable :: wtj      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r8_kind),    dimension(:,:,:), allocatable :: src_dist !< distance between destination grid and
                                                                        !! neighbor source grid.
   real(kind=r8_kind),    dimension(:,:), allocatable   :: rat_x    !< the ratio of coordinates of the dest grid
                                                                    !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                    !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r8_kind),    dimension(:,:), allocatable   :: rat_y  !< the ratio of coordinates of the dest grid
                                                                  !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                  !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r8_kind),    dimension(:), allocatable     :: lon_in   !< the coordinates of the source grid
   real(kind=r8_kind),    dimension(:), allocatable     :: lat_in   !< the coordinates of the source grid
   real(kind=r8_kind),    dimension(:), allocatable     :: area_frac_dst !< area fraction in destination grid.
   real(kind=r8_kind),    dimension(:,:), allocatable   :: mask_in
   real(kind=r8_kind)                                   :: max_src_dist
   logical                                              :: is_allocated !< set to true upon field allocation

end type horizInterpReals8_type

!> holds real(4) pointers for use in horiz_interp_type
type horizInterpReals4_type
   real(kind=r4_kind),    dimension(:,:), allocatable   :: faci     !< weights for conservative scheme
   real(kind=r4_kind),    dimension(:,:), allocatable   :: facj     !< weights for conservative scheme
   real(kind=r4_kind),    dimension(:,:), allocatable   :: area_src !< area of the source grid
   real(kind=r4_kind),    dimension(:,:), allocatable   :: area_dst !< area of the destination grid
   real(kind=r4_kind),    dimension(:,:,:), allocatable :: wti      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r4_kind),    dimension(:,:,:), allocatable :: wtj      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r4_kind),    dimension(:,:,:), allocatable :: src_dist !< distance between destination grid and
                                                                        !! neighbor source grid.
   real(kind=r4_kind),    dimension(:,:), allocatable   :: rat_x    !< the ratio of coordinates of the dest grid
                                                                    !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                    !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r4_kind),    dimension(:,:), allocatable   :: rat_y  !< the ratio of coordinates of the dest grid
                                                                  !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                  !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r4_kind),    dimension(:), allocatable     :: lon_in   !< the coordinates of the source grid
   real(kind=r4_kind),    dimension(:), allocatable     :: lat_in   !< the coordinates of the source grid
   real(kind=r4_kind),    dimension(:), allocatable     :: area_frac_dst !< area fraction in destination grid.
   real(kind=r4_kind),    dimension(:,:), allocatable   :: mask_in
   real(kind=r4_kind)                                   :: max_src_dist
   logical                                              :: is_allocated !< set to true upon field allocation

end type horizInterpReals4_type

!> Holds data pointers and metadata for horizontal interpolations, passed between the horiz_interp modules
!> @ingroup horiz_interp_type_mod
 type horiz_interp_type
   integer, dimension(:,:), allocatable   :: ilon    !< indices for conservative scheme
   integer, dimension(:,:), allocatable   :: jlat    !< indices for conservative scheme
                                                           !! wti ist used for derivative "weights" in bicubic
   integer, dimension(:,:,:), allocatable :: i_lon  !< indices for bilinear interpolation
                                                        !! and spherical regrid
   integer, dimension(:,:,:), allocatable :: j_lat  !< indices for bilinear interpolation
                                                        !! and spherical regrid
   logical, dimension(:,:), allocatable :: found_neighbors   !< indicate whether destination grid
                                                            !! has some source grid around it.
   integer, dimension(:,:), allocatable :: num_found
   integer                            :: nlon_src !< size of source grid
   integer                            :: nlat_src !< size of source grid
   integer                            :: nlon_dst !< size of destination grid
   integer                            :: nlat_dst !< size of destination grid
   integer                            :: interp_method      !< interpolation method.
                                                            !! =1, conservative scheme
                                                            !! =2, bilinear interpolation
                                                            !! =3, spherical regrid
                                                            !! =4, bicubic regrid
   logical                            :: I_am_initialized=.false.
   integer                            :: version                            !< indicate conservative
                                                                            !! interpolation version with value 1 or 2
   !--- The following are for conservative interpolation scheme version 2 ( through xgrid)
   integer                            :: nxgrid                             !< number of exchange grid
                                                                            !! between src and dst grid.
   integer, dimension(:), allocatable     :: i_src       !< indices in source grid.
   integer, dimension(:), allocatable     :: j_src       !< indices in source grid.
   integer, dimension(:), allocatable     :: i_dst       !< indices in destination grid.
   integer, dimension(:), allocatable     :: j_dst       !< indices in destination grid.
   type(horizInterpReals8_type) :: horizInterpReals8_type !< derived type holding kind 8 real data pointers
                                                                    !! if compiled with r8_kind
   type(horizInterpReals4_type) :: horizInterpReals4_type !< derived type holding kind 4 real data pointers
                                                                    !! if compiled with r8_kind
 end type

!> @addtogroup horiz_interp_type_mod
!> @{
contains

!######################################################################################################################
!> @brief horiz_interp_type_eq creates a copy of the horiz_interp_type object
 subroutine horiz_interp_type_eq(horiz_interp_out, horiz_interp_in)
    type(horiz_interp_type), intent(inout) :: horiz_interp_out !< Output object being set
    type(horiz_interp_type), intent(in)    :: horiz_interp_in !< Input object being copied

    if(.not.horiz_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'horiz_interp_type_eq: horiz_interp_type variable on right hand side is unassigned')
    endif

    horiz_interp_out%ilon            = horiz_interp_in%ilon
    horiz_interp_out%jlat            = horiz_interp_in%jlat
    horiz_interp_out%i_lon           = horiz_interp_in%i_lon
    horiz_interp_out%j_lat           = horiz_interp_in%j_lat
    horiz_interp_out%found_neighbors = horiz_interp_in%found_neighbors
    horiz_interp_out%num_found       = horiz_interp_in%num_found
    horiz_interp_out%nlon_src        =  horiz_interp_in%nlon_src
    horiz_interp_out%nlat_src        =  horiz_interp_in%nlat_src
    horiz_interp_out%nlon_dst        =  horiz_interp_in%nlon_dst
    horiz_interp_out%nlat_dst        =  horiz_interp_in%nlat_dst
    horiz_interp_out%interp_method   =  horiz_interp_in%interp_method
    horiz_interp_out%I_am_initialized = .true.
    horiz_interp_out%i_src           = horiz_interp_in%i_src
    horiz_interp_out%j_src           = horiz_interp_in%j_src
    horiz_interp_out%i_dst           = horiz_interp_in%i_dst
    horiz_interp_out%j_dst           = horiz_interp_in%j_dst

    if(horiz_interp_in%horizInterpReals8_type%is_allocated) then
      horiz_interp_out%horizInterpReals8_type%faci            = horiz_interp_in%horizInterpReals8_type%faci
      horiz_interp_out%horizInterpReals8_type%facj            = horiz_interp_in%horizInterpReals8_type%facj
      horiz_interp_out%horizInterpReals8_type%area_src        = horiz_interp_in%horizInterpReals8_type%area_src
      horiz_interp_out%horizInterpReals8_type%area_dst        = horiz_interp_in%horizInterpReals8_type%area_dst
      horiz_interp_out%horizInterpReals8_type%wti             = horiz_interp_in%horizInterpReals8_type%wti
      horiz_interp_out%horizInterpReals8_type%wtj             = horiz_interp_in%horizInterpReals8_type%wtj
      horiz_interp_out%horizInterpReals8_type%src_dist        = horiz_interp_in%horizInterpReals8_type%src_dist
      horiz_interp_out%horizInterpReals8_type%rat_x           = horiz_interp_in%horizInterpReals8_type%rat_x
      horiz_interp_out%horizInterpReals8_type%rat_y           = horiz_interp_in%horizInterpReals8_type%rat_y
      horiz_interp_out%horizInterpReals8_type%lon_in          = horiz_interp_in%horizInterpReals8_type%lon_in
      horiz_interp_out%horizInterpReals8_type%lat_in          = horiz_interp_in%horizInterpReals8_type%lat_in
      horiz_interp_out%horizInterpReals8_type%area_frac_dst   = horiz_interp_in%horizInterpReals8_type%area_frac_dst
      horiz_interp_out%horizInterpReals8_type%max_src_dist    =  horiz_interp_in%horizInterpReals8_type%max_src_dist
      horiz_interp_out%horizInterpReals8_type%is_allocated    = .true.
      ! this was left out previous to mixed mode
      horiz_interp_out%horizInterpReals8_type%mask_in         = horiz_interp_in%horizInterpReals8_type%mask_in

    else if (horiz_interp_in%horizInterpReals4_type%is_allocated) then
      horiz_interp_out%horizInterpReals4_type%faci            = horiz_interp_in%horizInterpReals4_type%faci
      horiz_interp_out%horizInterpReals4_type%facj            = horiz_interp_in%horizInterpReals4_type%facj
      horiz_interp_out%horizInterpReals4_type%area_src        = horiz_interp_in%horizInterpReals4_type%area_src
      horiz_interp_out%horizInterpReals4_type%area_dst        = horiz_interp_in%horizInterpReals4_type%area_dst
      horiz_interp_out%horizInterpReals4_type%wti             = horiz_interp_in%horizInterpReals4_type%wti
      horiz_interp_out%horizInterpReals4_type%wtj             = horiz_interp_in%horizInterpReals4_type%wtj
      horiz_interp_out%horizInterpReals4_type%src_dist        = horiz_interp_in%horizInterpReals4_type%src_dist
      horiz_interp_out%horizInterpReals4_type%rat_x           = horiz_interp_in%horizInterpReals4_type%rat_x
      horiz_interp_out%horizInterpReals4_type%rat_y           = horiz_interp_in%horizInterpReals4_type%rat_y
      horiz_interp_out%horizInterpReals4_type%lon_in          = horiz_interp_in%horizInterpReals4_type%lon_in
      horiz_interp_out%horizInterpReals4_type%lat_in          = horiz_interp_in%horizInterpReals4_type%lat_in
      horiz_interp_out%horizInterpReals4_type%area_frac_dst   = horiz_interp_in%horizInterpReals4_type%area_frac_dst
      horiz_interp_out%horizInterpReals4_type%max_src_dist    =  horiz_interp_in%horizInterpReals4_type%max_src_dist
      horiz_interp_out%horizInterpReals4_type%is_allocated    = .true.
      ! this was left out previous to mixed mode
      horiz_interp_out%horizInterpReals4_type%mask_in         = horiz_interp_in%horizInterpReals4_type%mask_in

    else
        call mpp_error(FATAL, "horiz_interp_type_eq: cannot assign unallocated real values from horiz_interp_in")
    endif

    if(horiz_interp_in%interp_method == CONSERVE) then
        horiz_interp_out%version =  horiz_interp_in%version
        if(horiz_interp_in%version==2) horiz_interp_out%nxgrid = horiz_interp_in%nxgrid
    end if

 end subroutine horiz_interp_type_eq
!######################################################################################################################

#include "horiz_interp_type_r4.fh"
#include "horiz_interp_type_r8.fh"

end module horiz_interp_type_mod
!> @}
! close documentation grouping
