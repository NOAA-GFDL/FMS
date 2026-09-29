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
!> @defgroup horiz_interp_type_mod horiz_interp_type_mod
!! @ingroup horiz_interp
!! @{
!! @author Zhi Liang
!!
!! @parblock
!! Defines horiz_interp_type
!! @endparblock

module horiz_interp_type_mod

use mpp_mod, only : mpp_send, mpp_recv, mpp_sync_self, mpp_error, FATAL
use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes
use mpp_mod, only : COMM_TAG_1, COMM_TAG_2
use platform_mod, only: r4_kind, r8_kind

implicit none
private


! parameter to determine interpolation method
 integer, parameter :: CONSERVE = 1 !< is an internally used parameter to mark conservative interpolation
 integer, parameter :: BILINEAR = 2 !< is an internally used parameter to mark bilinear interpolation
 integer, parameter :: SPHERICAL = 3 !< is an internally used parameter to mark spherical interpolation
 integer, parameter :: BICUBIC  = 4 !< is an internally used parameter to mark bicubic interpolation

public :: CONSERVE, BILINEAR, SPHERICAL, BICUBIC
public :: horiz_interp_type, stats, assignment(=)

!> Interface to override the "=" operator to a horiz_interp_eq that will
!! copy the Interp instance of horiz_interp_type with Interp2 = Interp1
interface assignment(=)
  module procedure horiz_interp_type_eq
end interface

!> Generic interfaces to compute statistics for bilinear and spherical interpolation
!! Calls stats_r4 when the output data is in 32-bit precision.  Calls stats_r8 when
!! the output data is in 64-bit precision.
interface stats
  module procedure stats_r4
  module procedure stats_r8
end interface


!> Nested derived type in horiz_interp_type that holds 64-bit interpolation weights, and
!! metadata.  If 64-bit grid and data arrays are provided to horiz_interp, Interp%horizInterpReals8_type
!! will be initialized and populated with mapping weights.
type horizInterpReals8_type
   real(kind=r8_kind), dimension(:,:), allocatable   :: faci
     !< holds weights for conservative interpolation version 1
   real(kind=r8_kind), dimension(:,:), allocatable   :: facj
     !< holds weights for conservative interpolation version 1
   real(kind=r8_kind), dimension(:,:), allocatable   :: area_src
     !< holds source grid area
   real(kind=r8_kind), dimension(:,:), allocatable   :: area_dst
     !< holds destination grid area.
   real(kind=r8_kind), dimension(:,:,:), allocatable :: wti
     !< holds interpolation weights for bilinear interpolation; x-derivatives for bicubic interpolation.
   real(kind=r8_kind), dimension(:,:,:), allocatable :: wtj
     !< holds interpolation weights for bilinear interpolation; y-derivatives for bicubic interpolation.
   real(kind=r8_kind), dimension(:,:,:), allocatable :: src_dist
     !< holds distance between destination grid and neighbor source grid in spherical interpolation.
   real(kind=r8_kind), dimension(:,:), allocatable   :: rat_x
     !< holds (x_dest -x_src_r)/(x_src_l -x_src_r) for bicubic interpolation
   real(kind=r8_kind), dimension(:,:), allocatable   :: rat_y
     !< holds (y_src_l -y_src_r)/(y_src_l -y_src_r) for bicubic interpolation
   real(kind=r8_kind), dimension(:), allocatable     :: lon_in
     !< holds the longitude coordinates on the source grid
   real(kind=r8_kind), dimension(:), allocatable     :: lat_in
     !< holds the latitude coordinates on the source grid
   real(kind=r8_kind), dimension(:), allocatable     :: area_frac_dst
     !< holds interpolation weights for conservative interpolation, version2
   real(kind=r8_kind), dimension(:,:), allocatable   :: mask_in
     !< masks the input grid to skip input cells when interpolation
   real(kind=r8_kind) :: max_src_dist
     !< sets to max_dist in spherical interpolation
   logical :: is_allocated = .false.
     !< is .true. if Interp is populated

end type horizInterpReals8_type

!> Nested derived type in horiz_interp_type that holds 32-bit interpolation weights, and
!! metadata.  If 32-bit grid and data arrays are provided to horiz_interp, Interp%horizInterpReals4_type
!! will be initialized and populated with mapping weights.
type horizInterpReals4_type
   real(kind=r4_kind), dimension(:,:), allocatable   :: faci
     !< holds weights for conservative interpolation version 1
   real(kind=r4_kind), dimension(:,:), allocatable   :: facj
     !< holds weights for conservative interpolation version 1
   real(kind=r4_kind), dimension(:,:), allocatable   :: area_src
     !< holds source grid area.
   real(kind=r4_kind), dimension(:,:), allocatable   :: area_dst
     !< holds destination grid area.
   real(kind=r4_kind), dimension(:,:,:), allocatable :: wti
     !< holds interpolation weights for bilinear interpolation; x-derivatives for bicubic interpolation.
   real(kind=r4_kind), dimension(:,:,:), allocatable :: wtj
     !< holds interpolation weights for bilinear interpolation; y-derivatives for bicubic interpolation.
   real(kind=r4_kind), dimension(:,:,:), allocatable :: src_dist
     !< holds distance between destination grid and neighbor source grid in spherical interpolation.
   real(kind=r4_kind), dimension(:,:), allocatable   :: rat_x
     !< holds (x_dest -x_src_r)/(x_src_l -x_src_r) for bicubic interpolation
   real(kind=r4_kind), dimension(:,:), allocatable   :: rat_y
     !< holds (y_src_l -y_src_r)/(y_src_l -y_src_r) for bicubic interpolation
   real(kind=r4_kind), dimension(:), allocatable     :: lon_in
     !< holds the longitude coordinates on the source grid
   real(kind=r4_kind), dimension(:), allocatable     :: lat_in
     !< holds the latitude coordinates on the source grid
   real(kind=r4_kind), dimension(:), allocatable     :: area_frac_dst
     !< holds interpolation weights for conservative interpolation, version2
   real(kind=r4_kind), dimension(:,:), allocatable   :: mask_in
     !< masks the input grid to skip input cells when interpolation
   real(kind=r4_kind) :: max_src_dist
     !< sets to max_dist in spherical interpolation
   logical :: is_allocated = .false.
     !< is .true. if Interp is populated

end type horizInterpReals4_type

!> Datatype holding interpolation weights, mapping indices, and metadata for horizontal interpolation.
!! All real members are stored in horizInterpReals8_type if the grid and data are represented in 64-bit
!! floating point precision or horizInterpReals4_type if the grid and data are represented in 32-bit
!! floating point precision.  Only one type, horizInterpReals4_type or horizInterpReals8_type, is allocated
!! and used for interpolation.
 type horiz_interp_type
   integer, dimension(:,:), allocatable   :: ilon
     !< contains the source grid mapping index in x-direction for conservative interpolation, version 1
   integer, dimension(:,:), allocatable   :: jlat
     !< contains the source grid mapping indices in y-direction for conservative interpolation, version 1
   integer, dimension(:,:,:), allocatable :: i_lon
   !< contains the source grid mapping indices in x-direction for bicubic and bilinear interpolation
   integer, dimension(:,:,:), allocatable :: j_lat
     !< contains the source grid cell indices in y-direction for bicubic and bilinear interpolation
   logical, dimension(:,:), allocatable :: found_neighbors
     !< is not used
   integer, dimension(:,:), allocatable :: num_found
     !< stores the number of neighbors found in spherical interpolation
   integer :: nlon_src
     !< is the size of source grid in the x direction
   integer :: nlat_src
     !< is the size of source grid in the y direction
   integer :: nlon_dst
     !< is the size of destination grid in the x direction
   integer :: nlat_dst
     !< is the size of destination grid in the y direction
   integer :: interp_method
     !< is the interpolation method set to 1 for conservative; 2 for bilinear; 3 for spherical; 4 for bicubic.
   logical :: I_am_initialized=.false.
     !< is set to .true. in horiz_interp_new.  Horiz_interp_base will fail if I_am_initialized = .false.
   integer :: version
     !< is set to 1 in horiz_interp_conserve_1dx1d_r4/8. Else set to 2; only used for conservative interpolation.
   integer :: nxgrid
     !< is the number of exchange grid cells for conservative interpolation, version 2.
   integer, dimension(:), allocatable :: i_src
     !< are the source grid mapping indices in the x-direction for conservative interpolation, version 2.
   integer, dimension(:), allocatable :: j_src
     !< are the source grid mapping indices in the y-direction for conservative interpolation, version 2.
   integer, dimension(:), allocatable :: i_dst
     !< are the destination grid mapping indices in the x-direction for conservative interpolation, version 2.
   integer, dimension(:), allocatable :: j_dst
     !< are the destination grid mapping indices in the y-direction for conservative interpolation, version 2.
   type(horizInterpReals8_type) :: horizInterpReals8_type
     !< holds more 64-bit floating point data required for interpolation.
   type(horizInterpReals4_type) :: horizInterpReals4_type
     !< holds more 32-bit floating point data required for interpolation.
 end type

contains

!######################################################################################################################
  !> @parblock
  !! Subroutine invoked when calling the "=" operator to copy all members of input horiz_interp_type
  !! into another instance of horiz_interp_type.  Do not call subroutine directly.  Instead, for copying,
  !! use the "=" operator:  Interp2 = Interp1.
  !! @endparblock
 subroutine horiz_interp_type_eq(horiz_interp_out, horiz_interp_in)
   type(horiz_interp_type), intent(inout) :: horiz_interp_out
     !< will contain the copied horiz_interp_type
   type(horiz_interp_type), intent(in) :: horiz_interp_in
     !< is the horiz_interp_type to copy

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
