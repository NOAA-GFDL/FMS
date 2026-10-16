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
!> @defgroup horiz_interp_spherical_mod horiz_interp_spherical_mod
!! @ingroup horiz_interp
!! @{
!! @brief Performs spatial interpolation between grids using inverse-distance-weighted scheme.
!! @parblock
!! Horiz_interp_spherical_mod contains methods called from horiz_interp_mod to interpolate data
!! from rectangular/tripolar grid to rectangular/tripolar grid using the inverse-distance-weighted
!! interpolation.  Users are recommended to use the top-level horiz_interp_mod with "interp" set to
!! "spherical" for the inverse-distance-weighted interpolation.
!! @endparblock
module horiz_interp_spherical_mod

  use platform_mod,          only : r4_kind, r8_kind
  use mpp_mod,               only : mpp_error, FATAL, WARNING, stdout
  use mpp_mod,               only : mpp_root_pe, mpp_pe
  use mpp_mod,               only : input_nml_file
  use fms_mod,               only : write_version_number
  use fms_mod,               only : check_nml_error
  use constants_mod,         only : pi
  use horiz_interp_type_mod, only : horiz_interp_type, stats, SPHERICAL

  implicit none
  private

  !> Generic interface to interpolate data from source grid to target grid from the weights stored in Interp.
  !! Calls horiz_interp_spherical_r4 if input and output data are 2D arrays in 32-bit precision.
  !! Calls horiz_interp_spherical_r8 if input and output data are 2D arrays in 64-bit precision.
  interface horiz_interp_spherical
    module procedure horiz_interp_spherical_r4
    module procedure horiz_interp_spherical_r8
  end interface

  !> Generic interface to compute interpolation weights and mapping indices.
  !! Calls horiz_interp_spherical_new_r4 when input and output 2D grid arrays are in 32-bit precision.
  !! Calls horiz_interp_spherical_new_r8 when input and output 2D grid arrays are in 64-bit precision.
  interface horiz_interp_spherical_new
    module procedure horiz_interp_spherical_new_r4
    module procedure horiz_interp_spherical_new_r8
  end interface

  !> Unused interface
  interface horiz_interp_spherical_wght
    module procedure horiz_interp_spherical_wght_r4
    module procedure horiz_interp_spherical_wght_r8
  end interface

  public :: horiz_interp_spherical_new, horiz_interp_spherical, horiz_interp_spherical_del
  public :: horiz_interp_spherical_init, horiz_interp_spherical_wght


  !> Generic interface used internally when search_method = "full_search" to search
  !! for nearest neighbors. Calls full_search_r4 when the input and output grids are in 32-bit
  !! precision.  Calls full_search_r8 when the input and output grids are in 64-bit precision.
  interface full_search
    module procedure full_search_r4
    module procedure full_search_r8
  end interface

  !> Generic interface used internally when search_method = "radial_search" to search
  !! for nearest neighbors. Calls radial_search_r4 when the input and output grids are
  !! in 32-bit precision.  Calls radial_search_r8 when the input and output grids are in
  !! 64-bit precision.
  interface radial_search
    module procedure radial_search_r4
    module procedure radial_search_r8
  end interface

  !> Generic inteface used internally in full_search and radial_search.
  !! Calls spherical_distance_r4 when the lon and lat values are in 32-bit precision.
  !! Calls spherical_distance_r8 when the lon and lat values are in 64-bit precision.
  interface spherical_distance
    module procedure spherical_distance_r4
    module procedure spherical_distance_r8
  end interface

  integer, parameter :: max_neighbors = 400
    !< is the maximum number of neighboring input grid cells around output grid cell.
  real(R8_KIND), parameter :: max_dist_default = 0.1_r8_kind
    !< is the maximum angular distance in radians for input and output grid cells to neighbor.
  integer, parameter :: num_nbrs_default = 4
    !< is the minimum number of neighbors to compute for each output grid cell.
  real(R8_KIND), parameter :: large=1.e20_r8_kind !< is a large number
  real(R8_KIND),  parameter :: epsln=1.e-10_r8_kind !< is a small number

  integer            :: pe, root_pe


  character(len=32) :: search_method = "radial_search"
    !< is a namelist flag to set the nearest neighbor searching method to either "radial_search" (default)
    !! or "full_search".  When set to "radial_search", the search may not be as accurate for some cases.
    !! Normally the search will be ok if you chose suitable max_dist. When search_method is "full_search",
    !! it will be always accurate, but will be slower comparing to "radial_search".  Normally these two search
    !! algorithm will produce same results other than order of operation. "radial_search" are recommended to use.
    !! The purpose to add "full_search" is in case you think you interpolation results is not right,
    !! you have other option to verify.

  namelist /horiz_interp_spherical_nml/ search_method

#include<file_version.h>
  logical            :: module_is_initialized = .FALSE.

contains

  !#######################################################################

  !> @parblock
  !! Initializes horiz_interp_spherical_mod.  Called from horiz_interp_init in horiz_interp_mod
  !! @endparblock
  subroutine horiz_interp_spherical_init
    integer :: ierr, io


    if(module_is_initialized) return
    call write_version_number("horiz_interp_spherical_mod", version)
    read (input_nml_file, horiz_interp_spherical_nml, iostat=io)
    ierr = check_nml_error(io,'horiz_interp_spherical_nml')

     module_is_initialized = .true.

  end subroutine horiz_interp_spherical_init

  !#######################################################################

  !> @parblock
  !! Deallocates arrays holding spherical interpolation weights and mapping indices
  !! in Interp.  Resets %is_allocated to .false.  Called from horiz_interp_del in horiz_interp_mod.
  !! @endparblock
  subroutine horiz_interp_spherical_del( Interp )

    type (horiz_interp_type), intent(inout) :: Interp !< will be reset with deallocated memory

    if(Interp%horizInterpReals4_type%is_allocated) then
      if(allocated(Interp%horizInterpReals4_type%src_dist)) deallocate(Interp%horizInterpReals4_type%src_dist)
    else if (Interp%horizInterpReals8_type%is_allocated) then
      if(allocated(Interp%horizInterpReals8_type%src_dist)) deallocate(Interp%horizInterpReals8_type%src_dist)
    endif
    if(allocated(Interp%num_found)) deallocate(Interp%num_found)
    if(allocated(Interp%i_lon))     deallocate(Interp%i_lon)
    if(allocated(Interp%j_lat))     deallocate(Interp%j_lat)

    Interp%horizInterpReals4_type%is_allocated = .false.
    Interp%horizInterpReals8_type%is_allocated = .false.

  end subroutine horiz_interp_spherical_del

  !#######################################################################

#include "horiz_interp_spherical_r4.fh"
#include "horiz_interp_spherical_r8.fh"

end module horiz_interp_spherical_mod
!> @}
! close documentation grouping
