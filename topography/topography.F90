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
!> @defgroup topography_mod topography_mod
!> @ingroup topography
!> @brief Routines for creating land surface topography fields and land-water masks
!! for latitude-longitude grids.
!> @author Bruce Wyman
!!
!! This module generates realistic mountains and land-water masks
!! on a specified latitude-longitude grid by interpolating from the
!! 1/6 degree Navy mean topography and percent water data sets.
!! The fields that can be generated are mean and standard deviation
!! of topography within the specified grid boxes; and land-ocean (or
!! water) mask and land-ocean (or water) fractional area.
!!
!! The interpolation scheme conserves the area-weighted average
!! of the input data by using module horiz_interp.
!!
!! The interfaces get_gaussian_topog and gaussian_topog_init are documented
!! in @ref gaussian_topog_mod

module topography_mod

use gaussian_topog_mod, only: gaussian_topog_init, get_gaussian_topog
use   horiz_interp_mod, only: horiz_interp_type, horiz_interp_new, &
                              horiz_interp, horiz_interp_del

use            fms_mod, only: check_nml_error, stdlog,    &
                              mpp_pe, mpp_root_pe, write_version_number, &
                              error_mesg, FATAL, NOTE, &
                              mpp_error
!  required for fms2_io
use        fms2_io_mod, only: read_data, FmsNetcdfFile_t, file_exists, open_file
!-----------------------------------------------------------------------

use      constants_mod, only: PI
use            mpp_mod, only: input_nml_file
use platform_mod,       only: r4_kind, r8_kind, FMS_PATH_LEN

implicit none
private

public :: topography_init,                 &
          get_topog_mean, get_topog_stdev, &
          get_ocean_frac, get_ocean_mask,  &
          get_water_frac, get_water_mask,  &
          gaussian_topog_init, get_gaussian_topog

!> @brief Returns a "realistic" mean surface height field.
!!
!> Returns realistic mountains on a latitude-longtude grid.
!! The returned field is the mean topography for the given grid boxes.
!! Computed using a conserving area-weighted interpolation.
!! The current input data set is the 1/6 degree Navy mean topography.
!!
!! @param blon The longitude (in radians) at grid box boundaries.
!! @param blat The latitude (in radians) at grid box boundaries.
!! @param zmean The mean surface height (meters). The size of this
!!                        field must be size(blon)-1 by size(blat)-1.
!! @return A logical value of TRUE is returned if the surface height field
!! was successfully created. A value of FALSE may be returned if the
!! input topography data set was not readable.
!!
!! <br>Example usage:
!! @code{.F90} flag = get_topog_mean ( blon, blat, zmean )@endcode
!> @ingroup topography_mod

interface get_topog_mean
    module procedure get_topog_mean_1d_r4, get_topog_mean_1d_r8
    module procedure get_topog_mean_2d_r4, get_topog_mean_2d_r8
end interface get_topog_mean

!> Returns a standard deviation of higher resolution topography with
!! the given model grid boxes.
!!
!> Returns the standard deviation of the "finer" input topography data set,
!! currently the Navy 1/6 degree mean topography data, within the
!! boundaries of the given input grid.
!!
!! @param blon The longitude (in radians) at grid box boundaries.
!! @param blat The latitude (in radians) at grid box boundaries.
!! @param [out] stdev The standard deviation of surface height (in meters) within
!! given input model grid boxes.
!! The size of this field must be size(blon)-1 by size(blat)-1.
!!
!! @return A logical value of TRUE is returned if the output field was
!! successfully created. A value of FALSE may be returned if the
!! input topography data set was not readable.
!!
!! Example usage:
!! @code{.F90} flag = get_topog_stdev( blon, blat, stdev ) @code
!> @ingroup topography_mod
interface get_topog_stdev
    module procedure get_topog_stdev_1d_r4, get_topog_stdev_1d_r8
    module procedure get_topog_stdev_2d_r4, get_topog_stdev_2d_r8
end interface get_topog_stdev

!> @brief Returns fractional area covered by ocean in a grid box.
!! Returns fractional area covered by ocean in the given model grid boxes.
!!
!! @param blon The longitude (in radians) at grid box boundaries.
!! @param blat The latitude (in radians) at grid box boundaries.
!! @param ocean_frac The fractional amount (0 to 1) of ocean in a grid box.
!! The size of this field must be size(blon)-1 by size(blat)-1.
!!
!! @return A logical value of TRUE is returned if the output field
!! was successfully created. A value of FALSE may be returned
!! if the Navy 1/6 degree percent water data set was not readable.
!!
!! Example usage:
!! @code{.F90} flag = get_ocean_frac ( blon, blat, ocean_frac ) @endcode
!> @ingroup topography_mod
interface get_ocean_frac
    module procedure get_ocean_frac_1d_r4, get_ocean_frac_1d_r8
    module procedure get_ocean_frac_2d_r4, get_ocean_frac_2d_r8
end interface get_ocean_frac

!> @brief Returns a land-ocean mask in a grid box.
!!
!> Returns a land-ocean mask in the given model grid boxes.
!!
!! @param blon The longitude (in radians) at grid box boundaries.
!! @param blat The latitude (in radians) at grid box boundaries.
!! @param ocean_frac The fractional amount (0 to 1) of ocean in a grid box.
!!  The size of this field must be size(blon)-1 by size(blat)-1.
!!
!! @return A logical value of TRUE is returned if the output field
!! was successfully created. A value of FALSE may be returned
!! if the Navy 1/6 degree percent water data set was not readable.
!!
!! Example code:
!! @code{.F90} flag = get_ocean_mask( blon, blat, ocean_mask ) @endcode
!> @ingroup topography_mod
interface get_ocean_mask
    module procedure get_ocean_mask_1d_r4, get_ocean_mask_1d_r8
    module procedure get_ocean_mask_2d_r4, get_ocean_mask_2d_r8
end interface get_ocean_mask

!> @brief Returns fractional area covered by water.
!!
!> Returns the percent of water in a grid box.
!!
!! @param blon The longitude (in radians) at grid box boundaries.
!! @param blat The latitude (in radians) at grid box boundaries.
!! @param [out] water_frac The fractional amount (0 to 1) of water in a grid box.
!!     The size of this field must be size(blon)-1 by size(blat)-1.
!!
!! @return A logical value of TRUE is returned if the output field
!! was successfully created. A value of FALSE may be returned
!! if the Navy 1/6 degree percent water data set was not readable.
!!
!! <br>Example usage:<br> @code{.F90} flag = get_water_frac ( blon, blat, water_frac ) @endcode
!> @ingroup topography_mod
interface get_water_frac
    module procedure get_water_frac_1d_r4, get_water_frac_1d_r8
    module procedure get_water_frac_2d_r4, get_water_frac_2d_r8
end interface get_water_frac

!> @brief Returns a land-water mask in a grid box.
!!
!> Returns a land-water mask in the given model grid boxes.
!!
!! @param blon The longitude (in radians) at grid box boundaries.
!! @param blat The latitude (in radians) at grid box boundaries.
!! @param water_mask A binary mask for water (true) or land (false).
!! The size of this field must be size(blon)-1 by size(blat)-1.
!!
!! @return A logical value of TRUE is returned if the output field
!! was successfully created. A value of FALSE may be returned
!! if the Navy 1/6 degree percent water data set was not readable.
!!
!! Example usage: @code{.F90}flag = get_water_mask( blon, blat, water_mask ) @endcode
!> @ingroup topography_mod
interface get_water_mask
    module procedure get_water_mask_1d_r4, get_water_mask_1d_r8
    module procedure get_water_mask_2d_r4, get_water_mask_2d_r8
end interface get_water_mask

interface interp_topog
    module procedure interp_topog_1d_r4, interp_topog_1d_r8
    module procedure interp_topog_2d_r4, interp_topog_2d_r8
end interface interp_topog

interface find_indices
    module procedure find_indices_r4, find_indices_r8
end interface find_indices

interface input_data
    module procedure input_data_r4, input_data_r8
end interface input_data

interface interp_water
    module procedure interp_water_1d_r4, interp_water_1d_r8
    module procedure interp_water_2d_r4, interp_water_2d_r8
end interface interp_water

interface determine_ocean_points
    module procedure determine_ocean_points_r4, determine_ocean_points_r8
end interface determine_ocean_points

!> @addtogroup topography_mod
!> @{

character(len=FMS_PATH_LEN) :: topog_file = 'DATA/navy_topography.data', &
                               water_file = 'DATA/navy_pctwater.data'
namelist /topography_nml/ topog_file, water_file

integer, parameter    :: TOPOG_INDEX = 1
integer, parameter    :: WATER_INDEX = 2
logical :: file_is_opened(2) = .false.
type(FmsNetcdfFile_t) :: fileobj(2) !< needed for fms2_io

!-----------------------------------------------------------------------
! --- resolution of the topography data set ---
! <DATASET NAME="">
!   This module uses the 1/6 degree U.S. Navy mean topography
!   and percent water data sets.
!
!   These data sets have been re-formatted to separate 32-bit IEEE files.
!   The names of these files is specified by the <LINK SRC="#NAMELIST">namelist</LINK> input.
!
!The format for both files is as follows:
! <PRE>
!     record = 1    nlon, nlat
!     record = 2    blon, blat
!     record = 3    data
! </PRE>
!where:
! <PRE>
!     nlon, nlat = The number of longitude and latitude points
!                  in the horizontal grid.  For the 1/6 degree
!                  data sets this is 2160 x 1080. [integer]
!     blon, blat = The longitude and latitude grid box boundaries in degrees.
!                     [real :: blon(nlon+1), blat(nlat+1)]
!
!     data       = The topography or percent water data.
!                    [real :: data(nlon,nlat)]
! </PRE>
! </DATASET>
integer :: ipts, jpts
integer, parameter :: COMPUTE_STDEV = 123  ! use this flag to
                                             !   compute st dev

!-----------------------------------------------------------------------

! Include variable "version" to be written to log file.
#include<file_version.h>

logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine topography_init ()
   if ( module_is_initialized ) return
   call write_version_number("TOPOGRAPHY_MOD", version)
   call read_namelist
   module_is_initialized = .TRUE.
end subroutine topography_init

!#######################################################################
!#######################################################################
!##################   private interfaces below here   ##################
!#######################################################################

function open_topog_file ( )
logical :: open_topog_file
real(kind=r4_kind)    :: r_ipts, r_jpts
integer :: namelen

namelen = len(trim(topog_file))
   if ( file_exists(topog_file) .AND. topog_file(namelen-2:namelen) == '.nc') then
      if (mpp_pe() == mpp_root_pe()) call mpp_error ('topography_mod', &
            'Reading NetCDF formatted input data file: '//trim(topog_file), NOTE)
         if(.not. file_is_opened(TOPOG_INDEX) ) then
            if(.not. open_file(fileobj(TOPOG_INDEX), topog_file, 'read' )) then
               call mpp_error(FATAL, 'topography_mod: Error in opening file '//trim(topog_file))
            endif
         endif

         call read_data(fileobj(TOPOG_INDEX), 'ipts', r_ipts)
         call read_data(fileobj(TOPOG_INDEX), 'jpts', r_jpts)
         ipts = nint(r_ipts)
         jpts = nint(r_jpts)
         open_topog_file = .true.
         file_is_opened(TOPOG_INDEX) = .true.
      else
         open_topog_file = .false.
      endif

end function open_topog_file

function open_water_file ( )
logical :: open_water_file
real(kind=r4_kind)    :: r_ipts, r_jpts
integer :: namelen

namelen = len(trim(water_file))
if ( file_exists(water_file) .AND. water_file(namelen-2:namelen) == '.nc') then
   if (mpp_pe() == mpp_root_pe()) call mpp_error ('topography_mod', &
      'Reading NetCDF formatted input data file: '//trim(water_file), NOTE)
      if(.not. file_is_opened(WATER_INDEX) ) then
         if(.not. open_file(fileobj(WATER_INDEX), water_file, 'read' )) then
            call mpp_error(FATAL, 'topography_mod: Error in opening file '//trim(water_file))
         endif
      endif

      call read_data(fileobj(WATER_INDEX), 'ipts', r_ipts)
      call read_data(fileobj(WATER_INDEX), 'jpts', r_jpts)
      ipts = nint(r_ipts)
      jpts = nint(r_jpts)
      open_water_file = .true.
      file_is_opened(WATER_INDEX) = .true.
   else
      open_water_file = .false.
   endif

end function open_water_file


!#######################################################################

!> @brief Reads the namelist file, write namelist to log file,
!! and initializes constants
subroutine read_namelist

integer :: iunit, ierr, io

!  read namelist

read (input_nml_file, topography_nml, iostat=io)
ierr = check_nml_error(io,'topography_nml')

!  write version and namelist to log file

if (mpp_pe() == mpp_root_pe()) then
   iunit = stdlog()
   write (iunit, nml=topography_nml)
endif

end subroutine read_namelist

#include "topography_r4.fh"
#include "topography_r8.fh"

end module topography_mod

! <INFO>

!   <TESTPROGRAM NAME="">
!
!  To run this program you will need the topography and percent water
!  data sets and use the following namelist (in the input nml file).
!
!   &amp;gaussian_topog_nml
!     height = 5000., 3000., 3000., 3000.,
!     olon   =   90.,  255.,  285.,    0.,
!     olat   =   45.,   45.,  -15.,  -90.,
!     wlon   =   15.,   10.,    5.,  180.,
!     wlat   =   15.,   25.,   25.,   20., /
!
!  program test
!
!  ! test program for topography and gaussian_topog modules
!  <PRE>
!  use topography_mod
!  implicit none
!
!  integer, parameter :: nlon=24, nlat=18
!  real :: x(nlon), y(nlat), xb(nlon+1), yb(nlat+1), z(nlon,nlat)
!  real :: hpi, rtd
!  integer :: i,j
!  logical :: a
!
!  ! gaussian mountain parameters
!  real, parameter :: ht=4000.
!  real, parameter :: x0=90., y0=45. ! origin in degrees
!  real, parameter :: xw=15., yw=15. ! half-width in degees
!  real, parameter :: xr=30., yr= 0. ! ridge-width in degrees
!
!  ! create lat/lon grid in radians
!    hpi = acos(0.0)
!    rtd = 90./hpi ! rad to deg
!    do i=1,nlon
!      xb(i) = 4.*hpi*real(i-1)/real(nlon)
!    enddo
!      xb(nlon+1) = xb(1)+4.*hpi
!      yb(1) = -hpi
!    do j=2,nlat
!      yb(j) = yb(j-1) + 2.*hpi/real(nlat)
!    enddo
!      yb(nlat+1) = hpi
!  ! mid-point of grid boxes
!    x(1:nlon) = 0.5*(xb(1:nlon)+xb(2:nlon+1))
!    y(1:nlat) = 0.5*(yb(1:nlat)+yb(2:nlat+1))
!  ! test topography_mod routines
!    a = get_topog_mean(xb,yb,z)
!    call printz ('get_topog_mean')
!
!    a = get_water_frac(xb,yb,z)
!    z = z*100. ! in percent
!    call printz ('get_water_frac')
!
!    a = get_ocean_frac(xb,yb,z)
!    z = z*100. ! in percent
!    call printz ('get_ocean_frac')
!
!  ! test gaussian_topog_mod routines
!    a = .true.
!    z = get_gaussian_topog(x,y,ht,x0,y0,xw,yw,xr,yr)
!    call printz ('get_gaussian_topog')
!
!    call gaussian_topog_init (x,y,z)
!    call printz ('gaussian_topog_init')
!
!  contains
!
!  ! simple printout of topog/water array
!    subroutine printz (lab)
!    character(len=*), intent(in) :: lab
!     if (a) then
!        print '(/a)', trim(lab)
!     else
!        print '(/a)', 'no data available: '//trim(lab)
!        return
!     endif
!      ! print full grid
!        print '(3x,25i5)', (nint(x(i)*rtd),i=1,nlon)
!      do j=nlat,1,-1
!        print '(i3,25i5)', nint(y(j)*rtd), (nint(z(i,j)),i=1,nlon)
!      enddo
!    end subroutine printz
!
!  end program test
!   </PRE>
!   </TESTPROGRAM>

!   <BUG>
!      Water mask produces some possible erroneous water points along
!      the coast of Antarctic (at about 90W).
!   </BUG>

!   <FUTURE>Use of netcdf data sets. </FUTURE>
!   <FUTURE>Incorporate other topography and ocean data sets. </FUTURE>
!
! </INFO>
!> @}
! close documentation grouping
