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
!> @defgroup gaussian_topog_mod gaussian_topog_mod
!> @ingroup topography
!> @brief Routines for creating Gaussian-shaped land surface topography
!! for latitude-longitude grids.
!> @author Bruce Wyman
!!
!! Interfaces generate simple Gaussian-shaped mountains from
!! parameters specified by either argument list or namelist input.
!! The mountain shapes are controlled by the height, half-width,
!! and ridge-width parameters.

!> @addtogroup gaussian_topog_mod
!> @{
module gaussian_topog_mod

use  fms_mod, only: check_nml_error,                 &
                    stdlog, write_version_number,    &
                    mpp_pe, mpp_root_pe,             &
                    error_mesg, FATAL

use constants_mod, only: pi

use mpp_mod,       only: input_nml_file
use platform_mod,  only: r4_kind, r8_kind

implicit none
private

public :: gaussian_topog_init, get_gaussian_topog

interface gaussian_topog_init
   module procedure gaussian_topog_init_r4
   module procedure gaussian_topog_init_r8
end interface gaussian_topog_init

interface get_gaussian_topog
    module procedure get_gaussian_topog_r4, get_gaussian_topog_r8
end interface get_gaussian_topog

!! Namelist information for gaussian_topog_nml
!!
!!     The variables in this namelist are only used when routine
!!     <TT>gaussian_topog_init</TT> is called.  The namelist variables
!!     are dimensioned (by 10), so that multiple mountains can be generated.
!!
!!     Internal parameter mxmtns = 10. By default no mountains are generated.
!!    </DATA>
!!
!!     NAMELIST FOR GENERATING GAUSSIAN MOUNTAINS
!!
!!  * multiple mountains can be generated
!!  * the final mountains are the sum of all
!!
!!       height = height in meters
!!       olon, olat = longitude,latitude origin              (degrees)
!!       rlon, rlat = longitude,latitude half-width of ridge (degrees)
!!       wlon, wlat = longitude,latitude half-width of tail  (degrees)
!!
!!       Note: For the standard gaussian mountain
!!             set rlon = rlat = 0 .
!!
!! <PRE>
!!
!!       height -->   ___________________________
!!                   /                           \
!!                  /              |              \
!!    gaussian     /               |               \
!!      sides --> /                |                \
!!               /               olon                \
!!         _____/                olat                 \______
!!
!!              |    |             |
!!              |<-->|<----------->|
!!              |wlon|    rlon     |
!!               wlat     rlat
!!
   integer, parameter :: maxmts = 10

   real(kind=r8_kind), dimension(maxmts) :: height = 0.0_r8_kind !< height in meters of the gaussian mountiains
   real(kind=r8_kind), dimension(maxmts) ::  olon  = 0.0_r8_kind !< longitude of mountain origins (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  olat  = 0.0_r8_kind !< Latitude  of mountain origins (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  wlon  = 0.0_r8_kind !< Longitude of half-width mountain trails (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  wlat  = 0.0_r8_kind !< Latitude of half-width mountain trails (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  rlon  = 0.0_r8_kind !< Longitude of half-width mountain ridges (degrees)
                                                                !! for "standard" gaussian mountain set, rlon/rlat = 0
   real(kind=r8_kind), dimension(maxmts) ::  rlat  = 0.0_r8_kind !< Latitude of half-width mountain ridges (degrees)
                                                                !! for "standard" gaussian mountain set, rlon/rlat = 0

   namelist /gaussian_topog_nml/ height, olon, olat, wlon, wlat, rlon, rlat
! </NAMELIST>

!-----------------------------------------------------------------------

! Include variable "version" to be written to log file.
#include<file_version.h>

logical :: do_nml = .true.
logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine read_namelist

   integer :: iunit, ierr, io

!>  read namelist

   read (input_nml_file, gaussian_topog_nml, iostat=io)
   ierr = check_nml_error(io,'gaussian_topog_nml')

!>  write version and namelist to log file

   if (mpp_pe() == mpp_root_pe()) then
      iunit = stdlog()
      write (iunit, nml=gaussian_topog_nml)
   endif

   do_nml = .false.

end subroutine read_namelist

!#######################################################################

#include "gaussian_topog_r4.fh"
#include "gaussian_topog_r8.fh"

end module gaussian_topog_mod

!> @}
! close documentation grouping
