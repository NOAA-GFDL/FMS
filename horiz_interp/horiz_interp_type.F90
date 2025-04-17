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
 integer, parameter :: SPHERICAL = 3
 integer, parameter :: BICUBIC  = 4

public :: CONSERVE, BILINEAR, SPHERICAL, BICUBIC
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
   logical                                              :: is_allocated = .false. !< set to true upon field allocation
 contains
   procedure :: get_faci => get_faci_r8
   procedure :: get_facj => get_facj_r8
   procedure :: get_area_src => get_area_src_r8
   procedure :: get_area_dst => get_area_dst_r8
   procedure :: get_wti => get_wti_r8
   procedure :: get_wtj => get_wtj_r8
   procedure :: get_src_dist => get_src_dist_r8
   procedure :: get_rat_x => get_rat_x_r8
   procedure :: get_rat_y => get_rat_y_r8
   procedure :: get_lon_in => get_lon_in_r8
   procedure :: get_lat_in => get_lat_in_r8
   procedure :: get_area_frac_dst => get_area_frac_dst_r8
   procedure :: get_mask_in => get_mask_in_r8
   procedure :: get_max_src_dist => get_max_src_dist_r8
   procedure :: get_is_allocated => get_is_allocated_r8

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
   logical                                              :: is_allocated = .false. !< set to true upon field allocation
 contains
   procedure :: get_faci => get_faci_r4
   procedure :: get_facj => get_facj_r4
   procedure :: get_area_src => get_area_src_r4
   procedure :: get_area_dst => get_area_dst_r4
   procedure :: get_wti => get_wti_r4
   procedure :: get_wtj => get_wtj_r4
   procedure :: get_src_dist => get_src_dist_r4
   procedure :: get_rat_x => get_rat_x_r4
   procedure :: get_rat_y => get_rat_y_r4
   procedure :: get_lon_in => get_lon_in_r4
   procedure :: get_lat_in => get_lat_in_r4
   procedure :: get_area_frac_dst => get_area_frac_dst_r4
   procedure :: get_mask_in => get_mask_in_r4
   procedure :: get_max_src_dist => get_max_src_dist_r4
   procedure :: get_is_allocated => get_is_allocated_r4

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
   type(horizInterpReals8_type)           :: horizInterpReals8_type !< derived type holding kind 8 real data pointers
                                                                    !! if compiled with r8_kind
   type(horizInterpReals4_type)           :: horizInterpReals4_type !< derived type holding kind 4 real data pointers
                                                                    !! if compiled with r8_kind
  contains
   procedure :: get_ilon
   procedure :: get_jlat
   procedure :: get_i_lon
   procedure :: get_j_lat
   procedure :: get_found_neighbors
   procedure :: get_num_found
   procedure :: get_nlon_src
   procedure :: get_nlat_src
   procedure :: get_nlon_dst
   procedure :: get_nlat_dst
   procedure :: get_interp_method
   procedure :: get_I_am_initialized
   procedure :: get_version
   procedure :: get_nxgrid
   procedure :: get_i_src
   procedure :: get_j_src
   procedure :: get_i_dst
   procedure :: get_j_dst
   procedure :: get_horizInterpReals8_type
   procedure :: get_horizInterpReals4_type
 end type

!> @addtogroup horiz_interp_type_mod
!> @{
contains

!######################################################################################################################
 function get_faci_r8(this) result(faci)
    class(horizInterpReals8_type)                   :: this
    real(kind=r8_kind), dimension(:,:), allocatable :: faci

    if (allocated(this%faci)) then
     faci = this%faci
    endif

 end function get_faci_r8

 function get_faci_r4(this) result(faci)
    class(horizInterpReals4_type)                   :: this
    real(kind=r4_kind), dimension(:,:), allocatable :: faci

    if (allocated(this%faci)) then
     faci = this%faci
    endif

 end function get_faci_r4

 function get_facj_r8(this) result(facj)
    class(horizInterpReals8_type)                   :: this
    real(kind=r8_kind), dimension(:,:), allocatable :: facj

    if (allocated(this%facj)) then
     facj = this%facj
    endif

 end function get_facj_r8

 function get_facj_r4(this) result(facj)
    class(horizInterpReals4_type)                   :: this
    real(kind=r4_kind), dimension(:,:), allocatable :: facj

    if (allocated(this%facj)) then
     facj = this%facj
    endif

 end function get_facj_r4

 function get_area_src_r8(this) result(area_src)
    class(horizInterpReals8_type)                   :: this
    real(kind=r8_kind), dimension(:,:), allocatable :: area_src

    if (allocated(this%area_src)) then
     area_src = this%area_src
    endif

 end function get_area_src_r8

 function get_area_src_r4(this) result(area_src)
    class(horizInterpReals4_type)                   :: this
    real(kind=r4_kind), dimension(:,:), allocatable :: area_src

    if (allocated(this%area_src)) then
     area_src = this%area_src
    endif

 end function get_area_src_r4

 function get_area_dst_r8(this) result(area_dst)
    class(horizInterpReals8_type)                   :: this
    real(kind=r8_kind), dimension(:,:), allocatable :: area_dst

    if (allocated(this%area_dst)) then
     area_dst = this%area_dst
    endif

 end function get_area_dst_r8

 function get_area_dst_r4(this) result(area_dst)
    class(horizInterpReals4_type)                   :: this
    real(kind=r4_kind), dimension(:,:), allocatable :: area_dst

    if (allocated(this%area_dst)) then
     area_dst = this%area_dst
    endif

 end function get_area_dst_r4

 function get_wti_r8(this) result(wti)
    class(horizInterpReals8_type)                     :: this
    real(kind=r8_kind), dimension(:,:,:), allocatable :: wti

    if (allocated(this%wti)) then
     wti = this%wti
    endif

 end function get_wti_r8

 function get_wti_r4(this) result(wti)
    class(horizInterpReals4_type)                     :: this
    real(kind=r4_kind), dimension(:,:,:), allocatable :: wti

    if (allocated(this%wti)) then
     wti = this%wti
    endif

 end function get_wti_r4

 function get_wtj_r8(this) result(wtj)
    class(horizInterpReals8_type)                     :: this
    real(kind=r8_kind), dimension(:,:,:), allocatable :: wtj

    if (allocated(this%wtj)) then
     wtj = this%wtj
    endif

 end function get_wtj_r8

 function get_wtj_r4(this) result(wtj)
    class(horizInterpReals4_type)                     :: this
    real(kind=r4_kind), dimension(:,:,:), allocatable :: wtj

    if (allocated(this%wtj)) then
     wtj = this%wtj
    endif

 end function get_wtj_r4

 function get_src_dist_r8(this) result(src_dist)
    class(horizInterpReals8_type)                     :: this
    real(kind=r8_kind), dimension(:,:,:), allocatable :: src_dist

    if (allocated(this%src_dist)) then
     src_dist = this%src_dist
    endif

 end function get_src_dist_r8

 function get_src_dist_r4(this) result(src_dist)
    class(horizInterpReals4_type)                     :: this
    real(kind=r4_kind), dimension(:,:,:), allocatable :: src_dist

    if (allocated(this%src_dist)) then
     src_dist = this%src_dist
    endif

 end function get_src_dist_r4

 function get_rat_x_r8(this) result(rat_x)
    class(horizInterpReals8_type)                   :: this
    real(kind=r8_kind), dimension(:,:), allocatable :: rat_x

    if (allocated(this%rat_x)) then
     rat_x = this%rat_x
    endif

 end function get_rat_x_r8

 function get_rat_x_r4(this) result(rat_x)
    class(horizInterpReals4_type)                   :: this
    real(kind=r4_kind), dimension(:,:), allocatable :: rat_x

    if (allocated(this%rat_x)) then
     rat_x = this%rat_x
    endif

 end function get_rat_x_r4

 function get_rat_y_r8(this) result(rat_y)
    class(horizInterpReals8_type)                   :: this
    real(kind=r8_kind), dimension(:,:), allocatable :: rat_y

    if (allocated(this%rat_y)) then
     rat_y = this%rat_y
    endif

 end function get_rat_y_r8

 function get_rat_y_r4(this) result(rat_y)
    class(horizInterpReals4_type)                   :: this
    real(kind=r4_kind), dimension(:,:), allocatable :: rat_y

    if (allocated(this%rat_y)) then
     rat_y = this%rat_y
    endif

 end function get_rat_y_r4

 function get_lon_in_r8(this) result(lon_in)
    class(horizInterpReals8_type)                 :: this
    real(kind=r8_kind), dimension(:), allocatable :: lon_in

    if (allocated(this%lon_in)) then
     lon_in = this%lon_in
    endif

 end function get_lon_in_r8

 function get_lon_in_r4(this) result(lon_in)
    class(horizInterpReals4_type)                 :: this
    real(kind=r4_kind), dimension(:), allocatable :: lon_in

    if (allocated(this%lon_in)) then
     lon_in = this%lon_in
    endif

 end function get_lon_in_r4

 function get_lat_in_r8(this) result(lat_in)
    class(horizInterpReals8_type)                 :: this
    real(kind=r8_kind), dimension(:), allocatable :: lat_in

    if (allocated(this%lat_in)) then
     lat_in = this%lat_in
    endif

 end function get_lat_in_r8

 function get_lat_in_r4(this) result(lat_in)
    class(horizInterpReals4_type)                 :: this
    real(kind=r4_kind), dimension(:), allocatable :: lat_in

    if (allocated(this%lat_in)) then
     lat_in = this%lat_in
    endif

 end function get_lat_in_r4

 function get_area_frac_dst_r8(this) result(area_frac_dst)
    class(horizInterpReals8_type)                 :: this
    real(kind=r8_kind), dimension(:), allocatable :: area_frac_dst

    if (allocated(this%area_frac_dst)) then
     area_frac_dst = this%area_frac_dst
    endif

 end function get_area_frac_dst_r8

 function get_area_frac_dst_r4(this) result(area_frac_dst)
    class(horizInterpReals4_type)                 :: this
    real(kind=r4_kind), dimension(:), allocatable :: area_frac_dst

    if (allocated(this%area_frac_dst)) then
     area_frac_dst = this%area_frac_dst
    endif

 end function get_area_frac_dst_r4

 function get_mask_in_r8(this) result(mask_in)
    class(horizInterpReals8_type)                   :: this
    real(kind=r8_kind), dimension(:,:), allocatable :: mask_in

    if (allocated(this%mask_in)) then
     mask_in = this%mask_in
    endif

 end function get_mask_in_r8

 function get_mask_in_r4(this) result(mask_in)
    class(horizInterpReals4_type)                   :: this
    real(kind=r4_kind), dimension(:,:), allocatable :: mask_in

    if (allocated(this%mask_in)) then
     mask_in = this%mask_in
    endif

 end function get_mask_in_r4

 function get_max_src_dist_r8(this) result(max_src_dist)
    class(horizInterpReals8_type) :: this
    real(kind=r8_kind)            :: max_src_dist

      max_src_dist = this%max_src_dist

 end function get_max_src_dist_r8

 function get_max_src_dist_r4(this) result(max_src_dist)
    class(horizInterpReals4_type) :: this
    real(kind=r4_kind)            :: max_src_dist

      max_src_dist = this%max_src_dist

 end function get_max_src_dist_r4

 function get_is_allocated_r8(this) result(is_allocated)
    class(horizInterpReals8_type) :: this
    logical                       :: is_allocated

      is_allocated = this%is_allocated

 end function get_is_allocated_r8

 function get_is_allocated_r4(this) result(is_allocated)
    class(horizInterpReals4_type) :: this
    logical                       :: is_allocated

      is_allocated = this%is_allocated

 end function get_is_allocated_r4

 function get_ilon(this) result(ilon)
    class(horiz_interp_type)             :: this
    integer, dimension(:,:), allocatable :: ilon

    if (allocated(this%ilon)) then
      ilon = this%ilon
    endif

 end function get_ilon

 function get_jlat(this) result(jlat)
    class(horiz_interp_type)             :: this
    integer, dimension(:,:), allocatable :: jlat

    if (allocated(this%jlat)) then
      jlat = this%jlat
    endif

 end function get_jlat

 function get_i_lon(this) result(i_lon)
    class(horiz_interp_type)               :: this
    integer, dimension(:,:,:), allocatable :: i_lon

    if (allocated(this%i_lon)) then
      i_lon = this%i_lon
    endif

 end function get_i_lon

 function get_j_lat(this) result(j_lat)
    class(horiz_interp_type)               :: this
    integer, dimension(:,:,:), allocatable :: j_lat

    if (allocated(this%j_lat)) then
      j_lat = this%j_lat
    endif

 end function get_j_lat

 function get_found_neighbors(this) result(found_neighbors)
    class(horiz_interp_type)             :: this
    logical, dimension(:,:), allocatable :: found_neighbors

    if (allocated(this%found_neighbors)) then
      found_neighbors = this%found_neighbors
    endif

 end function get_found_neighbors

 function get_num_found(this) result(num_found)
    class(horiz_interp_type)             :: this
    integer, dimension(:,:), allocatable :: num_found

    if (allocated(this%num_found)) then
      num_found = this%num_found
    endif

 end function get_num_found

 function get_nlon_src(this) result(nlon_src)
    class(horiz_interp_type) :: this
    integer                  :: nlon_src

    nlon_src = this%nlon_src

 end function get_nlon_src

 function get_nlat_src(this) result(nlat_src)
    class(horiz_interp_type) :: this
    integer                  :: nlat_src

    nlat_src = this%nlat_src

 end function get_nlat_src

 function get_nlon_dst(this) result(nlon_dst)
    class(horiz_interp_type) :: this
    integer                  :: nlon_dst

    nlon_dst = this%nlon_dst

 end function get_nlon_dst

 function get_nlat_dst(this) result(nlat_dst)
    class(horiz_interp_type) :: this
    integer                  :: nlat_dst

    nlat_dst = this%nlat_dst

 end function get_nlat_dst

 function get_interp_method(this) result(interp_method)
    class(horiz_interp_type) :: this
    integer                  :: interp_method

    interp_method = this%interp_method

 end function get_interp_method

 function get_I_am_initialized(this) result(I_am_initialized)
    class(horiz_interp_type) :: this
    logical                  :: I_am_initialized

    I_am_initialized = this%I_am_initialized

 end function get_I_am_initialized

 function get_version(this) result(version)
    class(horiz_interp_type) :: this
    integer                  :: version

    version = this%version

 end function get_version

 function get_nxgrid(this) result(nxgrid)
    class(horiz_interp_type) :: this
    integer                  :: nxgrid

    nxgrid = this%nxgrid

 end function get_nxgrid

 function get_i_src(this) result(i_src)
    class(horiz_interp_type)           :: this
    integer, dimension(:), allocatable :: i_src

    if (allocated(this%i_src)) then
      i_src = this%i_src
    endif

 end function get_i_src

 function get_j_src(this) result(j_src)
    class(horiz_interp_type)           :: this
    integer, dimension(:), allocatable :: j_src

    if (allocated(this%j_src)) then
      j_src = this%j_src
    endif

 end function get_j_src

 function get_i_dst(this) result(i_dst)
    class(horiz_interp_type)           :: this
    integer, dimension(:), allocatable :: i_dst

    if (allocated(this%i_dst)) then
      i_dst = this%i_dst
    endif

 end function get_i_dst

 function get_j_dst(this) result(j_dst)
    class(horiz_interp_type)           :: this
    integer, dimension(:), allocatable :: j_dst

    if (allocated(this%j_dst)) then
      j_dst = this%j_dst
    endif

 end function get_j_dst

 function get_horizInterpReals8_type(this) result(horizInterpReals8_type)
    class(horiz_interp_type)     :: this
    type(horizInterpReals8_type) :: horizInterpReals8_type

    !horizInterpReals8_type => this%horizInterpReals8_type

 end function get_horizInterpReals8_type

 function get_horizInterpReals4_type(this) result(horizInterpReals4_type)
    class(horiz_interp_type)     :: this
    type(horizInterpReals4_type) :: horizInterpReals4_type

    !horizInterpReals4_type => this%horizInterpReals4_type

 end function get_horizInterpReals4_type

!> @brief horiz_interp_type_eq creates a copy of the horiz_interp_type object
 subroutine horiz_interp_type_eq(horiz_interp_out, horiz_interp_in)
    type(horiz_interp_type), intent(inout) :: horiz_interp_out !< Output object being set
    type(horiz_interp_type), intent(in)    :: horiz_interp_in !< Input object being copied

    if(.not.horiz_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'horiz_interp_type_eq: horiz_interp_type variable on right hand side is unassigned')
    endif

    if( allocated(horiz_interp_in%ilon )) &
      horiz_interp_out%ilon = horiz_interp_in%ilon

    if( allocated(horiz_interp_in%jlat )) &
      horiz_interp_out%jlat = horiz_interp_in%jlat

    if( allocated(horiz_interp_in%i_lon )) &
      horiz_interp_out%i_lon = horiz_interp_in%i_lon

    if( allocated(horiz_interp_in%j_lat )) &
      horiz_interp_out%j_lat = horiz_interp_in%j_lat

    if( allocated(horiz_interp_in%found_neighbors )) &
      horiz_interp_out%found_neighbors = horiz_interp_in%found_neighbors

    if( allocated(horiz_interp_in%num_found )) &
      horiz_interp_out%num_found = horiz_interp_in%num_found

    if( allocated(horiz_interp_in%i_src )) &
      horiz_interp_out%i_src = horiz_interp_in%i_src

    if( allocated(horiz_interp_in%j_src )) &
      horiz_interp_out%j_src = horiz_interp_in%j_src

    if( allocated(horiz_interp_in%i_dst )) &
      horiz_interp_out%i_dst = horiz_interp_in%i_dst

    if( allocated(horiz_interp_in%j_dst )) &
      horiz_interp_out%j_dst = horiz_interp_in%j_dst

    horiz_interp_out%nlon_src =  horiz_interp_in%nlon_src
    horiz_interp_out%nlat_src =  horiz_interp_in%nlat_src
    horiz_interp_out%nlon_dst =  horiz_interp_in%nlon_dst
    horiz_interp_out%nlat_dst =  horiz_interp_in%nlat_dst
    horiz_interp_out%interp_method   =  horiz_interp_in%interp_method
    horiz_interp_out%I_am_initialized = .true.

    if(horiz_interp_in%horizInterpReals8_type%is_allocated) then

      if( allocated(horiz_interp_in%horizInterpReals8_type%faci)) &
        horiz_interp_out%horizInterpReals8_type%faci = horiz_interp_in%horizInterpReals8_type%faci

      if( allocated(  horiz_interp_in%horizInterpReals8_type%facj)) &
        horiz_interp_out%horizInterpReals8_type%facj = horiz_interp_in%horizInterpReals8_type%facj

      if( allocated(  horiz_interp_in%horizInterpReals8_type%area_src)) &
        horiz_interp_out%horizInterpReals8_type%area_src = horiz_interp_in%horizInterpReals8_type%area_src

      if( allocated(  horiz_interp_in%horizInterpReals8_type%area_dst)) &
        horiz_interp_out%horizInterpReals8_type%area_dst = horiz_interp_in%horizInterpReals8_type%area_dst

      if( allocated(  horiz_interp_in%horizInterpReals8_type%wti)) &
        horiz_interp_out%horizInterpReals8_type%wti = horiz_interp_in%horizInterpReals8_type%wti

      if( allocated(  horiz_interp_in%horizInterpReals8_type%wtj)) &
        horiz_interp_out%horizInterpReals8_type%wtj = horiz_interp_in%horizInterpReals8_type%wtj

      if( allocated(  horiz_interp_in%horizInterpReals8_type%src_dist)) &
        horiz_interp_out%horizInterpReals8_type%src_dist = horiz_interp_in%horizInterpReals8_type%src_dist

      if( allocated(  horiz_interp_in%horizInterpReals8_type%rat_x)) &
        horiz_interp_out%horizInterpReals8_type%rat_x = horiz_interp_in%horizInterpReals8_type%rat_x

      if( allocated(  horiz_interp_in%horizInterpReals8_type%rat_y)) &
        horiz_interp_out%horizInterpReals8_type%rat_y = horiz_interp_in%horizInterpReals8_type%rat_y

      if( allocated(  horiz_interp_in%horizInterpReals8_type%lon_in)) &
        horiz_interp_out%horizInterpReals8_type%lon_in = horiz_interp_in%horizInterpReals8_type%lon_in

      if( allocated(  horiz_interp_in%horizInterpReals8_type%lat_in)) &
        horiz_interp_out%horizInterpReals8_type%lat_in = horiz_interp_in%horizInterpReals8_type%lat_in

      if( allocated(  horiz_interp_in%horizInterpReals8_type%area_frac_dst)) &
        horiz_interp_out%horizInterpReals8_type%area_frac_dst = horiz_interp_in%horizInterpReals8_type%area_frac_dst

      horiz_interp_out%horizInterpReals8_type%max_src_dist =  horiz_interp_in%horizInterpReals8_type%max_src_dist

      horiz_interp_out%horizInterpReals8_type%is_allocated = .true.
      ! this was left out previous to mixed mode
      if( allocated(horiz_interp_in%horizInterpReals8_type%mask_in)) &
        horiz_interp_out%horizInterpReals8_type%mask_in = horiz_interp_in%horizInterpReals8_type%mask_in

    else if (horiz_interp_in%horizInterpReals4_type%is_allocated) then
      if( allocated(horiz_interp_in%horizInterpReals4_type%faci)) &
        horiz_interp_out%horizInterpReals4_type%faci = horiz_interp_in%horizInterpReals4_type%faci

      if( allocated(  horiz_interp_in%horizInterpReals4_type%facj)) &
        horiz_interp_out%horizInterpReals4_type%facj = horiz_interp_in%horizInterpReals4_type%facj

      if( allocated(  horiz_interp_in%horizInterpReals4_type%area_src)) &
        horiz_interp_out%horizInterpReals4_type%area_src = horiz_interp_in%horizInterpReals4_type%area_src

      if( allocated(  horiz_interp_in%horizInterpReals4_type%area_dst)) &
        horiz_interp_out%horizInterpReals4_type%area_dst = horiz_interp_in%horizInterpReals4_type%area_dst

      if( allocated(  horiz_interp_in%horizInterpReals4_type%wti)) &
        horiz_interp_out%horizInterpReals4_type%wti = horiz_interp_in%horizInterpReals4_type%wti

      if( allocated(  horiz_interp_in%horizInterpReals4_type%wtj)) &
        horiz_interp_out%horizInterpReals4_type%wtj = horiz_interp_in%horizInterpReals4_type%wtj

      if( allocated(  horiz_interp_in%horizInterpReals4_type%src_dist)) &
        horiz_interp_out%horizInterpReals4_type%src_dist = horiz_interp_in%horizInterpReals4_type%src_dist

      if( allocated(  horiz_interp_in%horizInterpReals4_type%rat_x)) &
        horiz_interp_out%horizInterpReals4_type%rat_x = horiz_interp_in%horizInterpReals4_type%rat_x

      if( allocated(  horiz_interp_in%horizInterpReals4_type%rat_y)) &
        horiz_interp_out%horizInterpReals4_type%rat_y = horiz_interp_in%horizInterpReals4_type%rat_y

      if( allocated(  horiz_interp_in%horizInterpReals4_type%lon_in)) &
        horiz_interp_out%horizInterpReals4_type%lon_in = horiz_interp_in%horizInterpReals4_type%lon_in

      if( allocated(  horiz_interp_in%horizInterpReals4_type%lat_in)) &
        horiz_interp_out%horizInterpReals4_type%lat_in = horiz_interp_in%horizInterpReals4_type%lat_in

      if( allocated(  horiz_interp_in%horizInterpReals4_type%area_frac_dst)) &
        horiz_interp_out%horizInterpReals4_type%area_frac_dst  = horiz_interp_in%horizInterpReals4_type%area_frac_dst

      horiz_interp_out%horizInterpReals4_type%max_src_dist =  horiz_interp_in%horizInterpReals4_type%max_src_dist

      horiz_interp_out%horizInterpReals4_type%is_allocated = .true.
      ! this was left out previous to mixed mode
      if( allocated(horiz_interp_in%horizInterpReals4_type%mask_in)) &
        horiz_interp_out%horizInterpReals4_type%mask_in = horiz_interp_in%horizInterpReals4_type%mask_in

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
