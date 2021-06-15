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

!> @file
!> @brief File for @ref gaussian_topog_mod

!> @addtogroup gaussian_topog_mod
!> @{
module gaussian_topog_mod

use  fms_mod, only: check_nml_error,                 &
                    stdlog, write_version_number,    &
                    mpp_pe, mpp_root_pe,             &
                    error_mesg, FATAL

use constants_mod, only: pi

use mpp_mod,       only: input_nml_file

implicit none
private

public :: gaussian_topog_init, get_gaussian_topog

!-----------------------------------------------------------------------
! <NAMELIST NAME="gaussian_topog_nml">
!   <DATA NAME="height" UNITS="meter" TYPE="real" DIM="(mxmtns)" DEFAULT="0.">
!     Height in meters of the Gaussian mountains.
!    </DATA>
!   <DATA NAME="olon, olat" UNITS="degree" TYPE="real" DIM="(mxmtns)" DEFAULT="0.">
!     The longitude and latitude of mountain origins (in degrees).
!    </DATA>
!   <DATA NAME="wlon, wlat" UNITS="degree" TYPE="real" DIM="(mxmtns)" DEFAULT="0.">
!     The longitude and latitude half-width of mountain tails (in degrees).
!    </DATA>
!   <DATA NAME="rlon, rlat" UNITS="degree" TYPE="real" DIM="(mxmtns)" DEFAULT="0.">
!     The longitude and latitude half-width of mountain ridges (in degrees).  For a
!     "standard" Gaussian mountain set rlon=rlat=0.
!    </DATA>
!
!    <DATA NAME="NOTE">
!     The variables in this namelist are only used when routine
!     <TT>gaussian_topog_init</TT> is called.  The namelist variables
!     are dimensioned (by 10), so that multiple mountains can be generated.
!
!     Internal parameter mxmtns = 10. By default no mountains are generated.
!    </DATA>

   integer, parameter :: maxmts = 10

   real, dimension(maxmts) :: height = 0.
   real, dimension(maxmts) ::  olon  = 0.
   real, dimension(maxmts) ::  olat  = 0.
   real, dimension(maxmts) ::  wlon  = 0.
   real, dimension(maxmts) ::  wlat  = 0.
   real, dimension(maxmts) ::  rlon  = 0.
   real, dimension(maxmts) ::  rlat  = 0.

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

!> Returns a surface height field that consists
!! of the sum of one or more Gaussian-shaped mountains.
!!
!> Returns a land surface topography that consists of a "set" of
!! simple Gaussian-shaped mountains.  The height, position,
!! width, and elongation of the mountains can be controlled
!! by variables in the namelist.
subroutine gaussian_topog_init ( lon, lat, zsurf )

real, intent(in)  :: lon(:) !< The mean grid box longitude in radians
real, intent(in)  :: lat(:) !< The mean grid box latitude in radians
real, intent(out) :: zsurf(:,:) !< The surface height (meters). Size must be size(lon) by size(lat)

integer :: n

  if (.not.module_is_initialized) then
     call write_version_number("GAUSSIAN_TOPOG_MOD", version)
  endif

  if(any(shape(zsurf) /= (/size(lon(:)),size(lat(:))/))) then
    call error_mesg ('get_gaussian_topog in topography_mod', &
     'shape(zsurf) is not equal to (/size(lon),size(lat)/)', FATAL)
  endif

  if (do_nml) call read_namelist

! compute sum of all non-zero mountains
  zsurf(:,:) = 0.
  do n = 1, maxmts
    if ( height(n) == 0. ) cycle
    zsurf = zsurf + get_gaussian_topog ( lon, lat, height(n), &
                olon(n), olat(n), wlon(n), wlat(n), rlon(n), rlat(n))
  enddo
 module_is_initialized = .TRUE.

end subroutine gaussian_topog_init

!#######################################################################

!> @brief Returns a simple surface height field that consists of a single
!! Gaussian-shaped mountain.
!!
!> The height, position, width, and elongation of the mountain
!! is controlled by optional arguments.
!! @param real lon The mean grid box longitude in radians.
!! @param real lat The mean grid box latitude in radians.
!! @param real height Maximum surface height in meters.
!! @param real olond, olatd Position/origin of mountain in degrees longitude and latitude.
!! This is the location of the maximum height.
!! @param real wlond, wlatd Gaussian half-width of mountain in degrees longitude and latitude.
!! @param real rlond, rlatd Ridge half-width of mountain in degrees longitude and latitude.
!! This is the elongation of the maximum height.
!! @param real zsurf The surface height (in meters).
!! The size of the returned field is size(lon) by size(lat).
!!   </OUT>
!!
!! @throws FATAL shape(zsurf) is not equal to (/size(lon),size(lat)/)
!!     Check the input grid size and output field size.
!!     The input grid is defined at the midpoint of grid boxes.
!!
!! @note
!!     Mountains do not wrap around the poles.
!
!! <br>Example usage:
!! @code{.F90} zsurf = <B>get_gaussian_topog</B> ( lon, lat, height
!!                    [, olond, olatd, wlond, wlatd, rlond, rlatd ] )@endcode
function get_gaussian_topog ( lon, lat, height,                          &
                              olond, olatd, wlond, wlatd, rlond, rlatd ) &
                     result ( zsurf )

real, intent(in)  :: lon(:), lat(:)
real, intent(in)  :: height
real, intent(in), optional :: olond, olatd, wlond, wlatd, rlond, rlatd
real :: zsurf(size(lon,1),size(lat,1))

integer :: i, j
real    :: olon, olat, wlon, wlat, rlon, rlat
real    :: tpi, dtr, dx, dy, xx, yy

  if (do_nml) call read_namelist

! no need to compute mountain if height=0
  if ( height == 0. ) then
       zsurf(:,:) = 0.
       return
  endif

  tpi = 2.0*pi
  dtr = tpi/360.

! defaults and convert degrees to radians (dtr)
  olon = 90.*dtr;  if (present(olond)) olon=olond*dtr
  olat = 45.*dtr;  if (present(olatd)) olat=olatd*dtr
  wlon = 15.*dtr;  if (present(wlond)) wlon=wlond*dtr
  wlat = 15.*dtr;  if (present(wlatd)) wlat=wlatd*dtr
  rlon =  0.    ;  if (present(rlond)) rlon=rlond*dtr
  rlat =  0.    ;  if (present(rlatd)) rlat=rlatd*dtr

! compute gaussian-shaped mountain
    do j=1,size(lat(:))
      dy = abs(lat(j) - olat)   ! dist from y origin
      yy = max(0., dy-rlat)/wlat
      do i=1,size(lon(:))
        dx = abs(lon(i) - olon) ! dist from x origin
        dx = min(dx, abs(dx-tpi))  ! To ensure that: -pi <= dx <= pi
        xx = max(0., dx-rlon)/wlon
        zsurf(i,j) = height*exp(-xx**2 - yy**2)
      enddo
    enddo

end function get_gaussian_topog

!#######################################################################

subroutine read_namelist

   integer :: unit, ierr, io
   real    :: dtr

!>  read namelist

   read (input_nml_file, gaussian_topog_nml, iostat=io)
   ierr = check_nml_error(io,'gaussian_topog_nml')

!>  write version and namelist to log file

   if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog()
      write (unit, nml=gaussian_topog_nml)
   endif

   do_nml = .false.

end subroutine read_namelist

!#######################################################################

end module gaussian_topog_mod

! <INFO>
!   <NOTE>
!     NAMELIST FOR GENERATING GAUSSIAN MOUNTAINS
!
!  * multiple mountains can be generated
!  * the final mountains are the sum of all
!
!       height = height in meters
!       olon, olat = longitude,latitude origin              (degrees)
!       rlon, rlat = longitude,latitude half-width of ridge (degrees)
!       wlon, wlat = longitude,latitude half-width of tail  (degrees)
!
!       Note: For the standard gaussian mountain
!             set rlon = rlat = 0 .
!
! <PRE>
!
!       height -->   ___________________________
!                   /                           \
!                  /              |              \
!    gaussian     /               |               \
!      sides --> /                |                \
!               /               olon                \
!         _____/                olat                 \______
!
!              |    |             |
!              |<-->|<----------->|
!              |wlon|    rlon     |
!               wlat     rlat
!
! </PRE>
!
!See the <LINK SRC="topography.html#TEST PROGRAM">topography </LINK>module documentation for a test program.
!   </NOTE>
! </INFO>
!> @}
! close documentation grouping
