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
!> @defgroup horiz_interp_spherical_mod horiz_interp_spherical_mod
!> @ingroup horiz_interp
!> @brief Performs spatial interpolation between grids using inverse-distance-weighted scheme.
!> This module can interpolate data from rectangular/tripolar grid
!! to rectangular/tripolar grid. The interpolation scheme is inverse-distance-weighted
!! scheme.    There is an optional mask field for missing input data.
!! An optional output mask field may be used in conjunction with
!! the input mask to show where output data exists.

!> @addtogroup horiz_interp_spherical_mod
!> @{
module horiz_interp_spherical_mod

  use platform_mod,          only : r4_kind, r8_kind
  use mpp_mod,               only : mpp_error, FATAL, WARNING, stdout
  use mpp_mod,               only : mpp_root_pe, mpp_pe
  use mpp_mod,               only : input_nml_file
  use fms_mod,               only : write_version_number
  use fms_mod,               only : check_nml_error
  use constants_mod,         only : pi
  use horiz_interp_type_mod, only : horiz_interp_type, stats

  implicit none
  private

  interface horiz_interp_spherical
    module procedure horiz_interp_spherical_r4
    module procedure horiz_interp_spherical_r8
  end interface

  interface horiz_interp_spherical_new
    module procedure horiz_interp_spherical_new_r4
    module procedure horiz_interp_spherical_new_r8
  end interface

  interface horiz_interp_spherical_wght
    module procedure horiz_interp_spherical_wght_r4
    module procedure horiz_interp_spherical_wght_r8
  end interface

  public :: horiz_interp_spherical_new, horiz_interp_spherical, horiz_interp_spherical_del
  public :: horiz_interp_spherical_init, horiz_interp_spherical_wght

 !> private helper routines
  interface full_search
    module procedure full_search_r4
    module procedure full_search_r8
  end interface

  interface radial_search
    module procedure radial_search_r4
    module procedure radial_search_r8
  end interface

  interface spherical_distance
    module procedure spherical_distance_r4
    module procedure spherical_distance_r8
  end interface

  integer, parameter :: max_neighbors = 400
  real(R8_KIND),    parameter :: max_dist_default = 0.1_r8_kind  ! radians
  integer, parameter :: num_nbrs_default = 4
  real(R8_KIND),    parameter :: large=1.e20_r8_kind
  real(R8_KIND),    parameter :: epsln=1.e-10_r8_kind

  integer            :: pe, root_pe


  character(len=32) :: search_method = "radial_search" !< Namelist variable to indicate the searching
                    !! method to find the
                    !! nearest neighbor points. Its value can be "radial_search" and "full_search",
                    !! with default value "radial_search". when search_method is "radial_search",
                    !! the search may be not quite accurate for some cases. Normally the search will
                    !! be ok if you chose suitable max_dist. When search_method is "full_search",
                    !! it will be always accurate, but will be slower comparing to "radial_search".
                    !! Normally these two search algorithm will produce same results other than
                    !! order of operation. "radial_search" are recommended to use. The purpose to
                    !! add "full_search" is in case you think you interpolation results is
                    !! not right, you have other option to verify.

!or "full_search"
  namelist /horiz_interp_spherical_nml/ search_method

  !-----------------------------------------------------------------------
  ! Include variable "version" to be written to log file.
#include<file_version.h>
  logical            :: module_is_initialized = .FALSE.

contains

  !#######################################################################

  !> Initializes module and writes version number to logfile.out
  subroutine horiz_interp_spherical_init
    integer :: ierr, io


    if(module_is_initialized) return
    call write_version_number("horiz_interp_spherical_mod", version)
    read (input_nml_file, horiz_interp_spherical_nml, iostat=io)
    ierr = check_nml_error(io,'horiz_interp_spherical_nml')

     module_is_initialized = .true.

  end subroutine horiz_interp_spherical_init

  !#######################################################################

  !> Deallocates memory used by "HI_KIND_TYPE" variables.
  !! Must be called before reinitializing with horiz_interp_spherical_new.
  subroutine horiz_interp_spherical_del( Interp )

    type (horiz_interp_type), intent(inout) :: Interp !< A derived-type variable returned by previous
                                           !! call to horiz_interp_spherical_new. The input variable
                                           !! must have allocated arrays. The returned variable will
                                           !! contain deallocated arrays.

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
