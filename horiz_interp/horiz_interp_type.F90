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
use platform_mod

implicit none
private


! parameter to determine interpolation method
 integer, parameter :: CONSERVE = 1
 integer, parameter :: BILINEAR = 2
 integer, parameter :: SPHERICA = 3
 integer, parameter :: BICUBIC  = 4

public :: CONSERVE, BILINEAR, SPHERICA, BICUBIC
public :: horiz_interp_type, stats, assignment(=)
public :: horiz_interp_r8_type, horiz_interp_r4_type

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

!<PUBLICTYPE >
!> @ingroup horiz_interp_type_mod
 type, abstract :: horiz_interp_type
   integer, dimension(:,:), pointer   :: ilon =>NULL()   !< indices for conservative scheme
   integer, dimension(:,:), pointer   :: jlat =>NULL()   !< indices for conservative scheme
   integer, dimension(:,:,:), pointer :: i_lon =>NULL() !< indices for bilinear interpolation
                                                        !! and spherical regrid
   integer, dimension(:,:,:), pointer :: j_lat =>NULL() !< indices for bilinear interpolation
                                                        !! and spherical regrid
   logical, dimension(:,:), pointer   :: found_neighbors =>NULL()       !< indicate whether destination grid
                                                                        !! has some source grid around it.
   integer, dimension(:,:), pointer   :: num_found => NULL()
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
   integer, dimension(:), pointer     :: i_src=>NULL()       !< indices in source grid.
   integer, dimension(:), pointer     :: j_src=>NULL()       !< indices in source grid.
   integer, dimension(:), pointer     :: i_dst=>NULL()       !< indices in destination grid.
   integer, dimension(:), pointer     :: j_dst=>NULL()       !< indices in destination grid.
 end type
!</PUBLICTYPE>

type, extends(horiz_interp_type) :: horiz_interp_r8_type
   real(r8_kind),    dimension(:,:), pointer   :: faci =>NULL()   !< weights for conservative scheme
   real(r8_kind),    dimension(:,:), pointer   :: facj =>NULL()   !< weights for conservative scheme
   real(r8_kind),    dimension(:,:), pointer   :: area_src =>NULL()              !< area of the source grid
   real(r8_kind),    dimension(:,:), pointer   :: area_dst =>NULL()              !< area of the destination grid
   real(r8_kind),    dimension(:,:,:), pointer :: wti =>NULL()      !< weights for bilinear interpolation
                                                           !! wti ist used for derivative "weights" in bicubic
   real(r8_kind),    dimension(:,:,:), pointer :: wtj =>NULL()      !< weights for bilinear interpolation
                                                           !! wti ist used for derivative "weights" in bicubic
   real(r8_kind),    dimension(:,:,:), pointer :: src_dist =>NULL()              !< distance between destination grid and
                                                                        !! neighbor source grid.
   real(r8_kind),    dimension(:,:), pointer   :: rat_x =>NULL() !< the ratio of coordinates of the dest grid
                                                        !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                        !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(r8_kind),    dimension(:,:), pointer   :: rat_y =>NULL() !< the ratio of coordinates of the dest grid
                                                        !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                        !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(r8_kind),    dimension(:), pointer     :: lon_in =>NULL()  !< the coordinates of the source grid
   real(r8_kind),    dimension(:), pointer     :: lat_in =>NULL()  !< the coordinates of the source grid
   real(r8_kind),    dimension(:), pointer     :: area_frac_dst=>NULL()              !< area fraction in destination grid.
   real(r8_kind),    dimension(:,:), pointer   :: mask_in=>NULL()
   real(r8_kind)                               :: max_src_dist
end type

type, extends(horiz_interp_type) :: horiz_interp_r4_type
   real(r4_kind),    dimension(:,:), pointer   :: faci =>NULL()   !< weights for conservative scheme
   real(r4_kind),    dimension(:,:), pointer   :: facj =>NULL()   !< weights for conservative scheme
   real(r4_kind),    dimension(:,:), pointer   :: area_src =>NULL()              !< area of the source grid
   real(r4_kind),    dimension(:,:), pointer   :: area_dst =>NULL()              !< area of the destination grid
   real(r4_kind),    dimension(:,:,:), pointer :: wti =>NULL()      !< weights for bilinear interpolation
                                                           !! wti ist used for derivative "weights" in bicubic
   real(r4_kind),    dimension(:,:,:), pointer :: wtj =>NULL()      !< weights for bilinear interpolation
                                                           !! wti ist used for derivative "weights" in bicubic
   real(r4_kind),    dimension(:,:,:), pointer :: src_dist =>NULL()              !< distance between destination grid and
                                                                        !! neighbor source grid.
   real(r4_kind),    dimension(:,:), pointer   :: rat_x =>NULL() !< the ratio of coordinates of the dest grid
                                                        !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                        !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(r4_kind),    dimension(:,:), pointer   :: rat_y =>NULL() !< the ratio of coordinates of the dest grid
                                                        !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                        !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(r4_kind),    dimension(:), pointer     :: lon_in =>NULL()  !< the coordinates of the source grid
   real(r4_kind),    dimension(:), pointer     :: lat_in =>NULL()  !< the coordinates of the source grid
   real(r4_kind),    dimension(:), pointer     :: area_frac_dst=>NULL()              !< area fraction in destination grid.
   real(r4_kind),    dimension(:,:), pointer   :: mask_in=>NULL()
   real(r4_kind)                               :: max_src_dist
end type

!> @addtogroup horiz_interp_type_mod
!> @{
contains

!######################################################################################################################
 subroutine horiz_interp_type_eq(horiz_interp_out, horiz_interp_in)
    class(horiz_interp_type), intent(inout) :: horiz_interp_out
    class(horiz_interp_type), intent(in)    :: horiz_interp_in
    logical :: valid_types
    valid_types = .false.

    if(.not.horiz_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'horiz_interp_type_eq: horiz_interp_type variable on right hand side is unassigned')
    endif
   select type (horiz_interp_out)
   type is(horiz_interp_r4_type)
   select type(horiz_interp_in)
   type is (horiz_interp_r4_type)
     horiz_interp_out%faci            => horiz_interp_in%faci
     horiz_interp_out%facj            => horiz_interp_in%facj
     horiz_interp_out%ilon            => horiz_interp_in%ilon
     horiz_interp_out%jlat            => horiz_interp_in%jlat
     horiz_interp_out%area_src        => horiz_interp_in%area_src
     horiz_interp_out%area_dst        => horiz_interp_in%area_dst
     horiz_interp_out%wti             => horiz_interp_in%wti
     horiz_interp_out%wtj             => horiz_interp_in%wtj
     horiz_interp_out%i_lon           => horiz_interp_in%i_lon
     horiz_interp_out%j_lat           => horiz_interp_in%j_lat
     horiz_interp_out%src_dist        => horiz_interp_in%src_dist
     horiz_interp_out%found_neighbors => horiz_interp_in%found_neighbors
     horiz_interp_out%max_src_dist    =  horiz_interp_in%max_src_dist
     horiz_interp_out%num_found       => horiz_interp_in%num_found
     horiz_interp_out%nlon_src        =  horiz_interp_in%nlon_src
     horiz_interp_out%nlat_src        =  horiz_interp_in%nlat_src
     horiz_interp_out%nlon_dst        =  horiz_interp_in%nlon_dst
     horiz_interp_out%nlat_dst        =  horiz_interp_in%nlat_dst
     horiz_interp_out%interp_method   =  horiz_interp_in%interp_method
     horiz_interp_out%rat_x           => horiz_interp_in%rat_x
     horiz_interp_out%rat_y           => horiz_interp_in%rat_y
     horiz_interp_out%lon_in          => horiz_interp_in%lon_in
     horiz_interp_out%lat_in          => horiz_interp_in%lat_in
     horiz_interp_out%I_am_initialized = .true.
     horiz_interp_out%i_src           => horiz_interp_in%i_src
     horiz_interp_out%j_src           => horiz_interp_in%j_src
     horiz_interp_out%i_dst           => horiz_interp_in%i_dst
     horiz_interp_out%j_dst           => horiz_interp_in%j_dst
     horiz_interp_out%area_frac_dst   => horiz_interp_in%area_frac_dst
     valid_types = .true.
   end select
   type is(horiz_interp_r8_type)
   select type(horiz_interp_in)
   type is (horiz_interp_r8_type)
     horiz_interp_out%faci            => horiz_interp_in%faci
     horiz_interp_out%facj            => horiz_interp_in%facj
     horiz_interp_out%ilon            => horiz_interp_in%ilon
     horiz_interp_out%jlat            => horiz_interp_in%jlat
     horiz_interp_out%area_src        => horiz_interp_in%area_src
     horiz_interp_out%area_dst        => horiz_interp_in%area_dst
     horiz_interp_out%wti             => horiz_interp_in%wti
     horiz_interp_out%wtj             => horiz_interp_in%wtj
     horiz_interp_out%i_lon           => horiz_interp_in%i_lon
     horiz_interp_out%j_lat           => horiz_interp_in%j_lat
     horiz_interp_out%src_dist        => horiz_interp_in%src_dist
     horiz_interp_out%found_neighbors => horiz_interp_in%found_neighbors
     horiz_interp_out%max_src_dist    =  horiz_interp_in%max_src_dist
     horiz_interp_out%num_found       => horiz_interp_in%num_found
     horiz_interp_out%nlon_src        =  horiz_interp_in%nlon_src
     horiz_interp_out%nlat_src        =  horiz_interp_in%nlat_src
     horiz_interp_out%nlon_dst        =  horiz_interp_in%nlon_dst
     horiz_interp_out%nlat_dst        =  horiz_interp_in%nlat_dst
     horiz_interp_out%interp_method   =  horiz_interp_in%interp_method
     horiz_interp_out%rat_x           => horiz_interp_in%rat_x
     horiz_interp_out%rat_y           => horiz_interp_in%rat_y
     horiz_interp_out%lon_in          => horiz_interp_in%lon_in
     horiz_interp_out%lat_in          => horiz_interp_in%lat_in
     horiz_interp_out%I_am_initialized = .true.
     horiz_interp_out%i_src           => horiz_interp_in%i_src
     horiz_interp_out%j_src           => horiz_interp_in%j_src
     horiz_interp_out%i_dst           => horiz_interp_in%i_dst
     horiz_interp_out%j_dst           => horiz_interp_in%j_dst
     horiz_interp_out%area_frac_dst   => horiz_interp_in%area_frac_dst
     valid_types = .true.
   end select
   end select
   if( .not. valid_types) call mpp_error(FATAL, "horiz_interp_type: invalid r4/r8 types passed to stats")
    if(horiz_interp_in%interp_method == CONSERVE) then
       horiz_interp_out%version =  horiz_interp_in%version
       if(horiz_interp_in%version==2) horiz_interp_out%nxgrid = horiz_interp_in%nxgrid
    end if

 end subroutine horiz_interp_type_eq
!######################################################################################################################
#undef FMS_HI_KIND
#define FMS_HI_KIND 4
#undef HI_STATS
#define HI_STATS stats_r4
#include <horiz_interp_type.inc> 
#undef FMS_HI_KIND
#define FMS_HI_KIND 8
#undef HI_STATS
#define HI_STATS stats_r8
#include <horiz_interp_type.inc> 

end module horiz_interp_type_mod
!> @}
! close documentation grouping
