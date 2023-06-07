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

function GET_GAUSSIAN_TOPOG_(lon, lat, height,                         &
                            olond, olatd, wlond, wlatd, rlond, rlatd ) &
                     result (zsurf )

real(kind=FMS_TOP_KIND_), intent(in)  :: lon(:), lat(:)
real(kind=FMS_TOP_KIND_), intent(in)  :: height
real(kind=FMS_TOP_KIND_), intent(in), optional :: olond, olatd, wlond, wlatd, rlond, rlatd
real(kind=FMS_TOP_KIND_) :: zsurf(size(lon,1),size(lat,1))

integer :: i, j
real(kind=FMS_TOP_KIND_)    :: olon, olat, wlon, wlat, rlon, rlat
real(kind=FMS_TOP_KIND_)    :: tpi, dtr, dx, dy, xx, yy
integer, parameter     :: lkind = FMS_TOP_KIND_

  if (do_nml) call read_namelist

! no need to compute mountain if height=0
  if ( height == 0.0_lkind) then
       zsurf(:,:) = 0.0_lkind
       return
  endif

  tpi = 2.0_lkind*real(pi, FMS_TOP_KIND_)
  dtr = tpi/360.0_lkind

! defaults and convert degrees to radians (dtr)
  olon = 90.0_lkind*dtr;  if (present(olond)) olon=olond*dtr
  olat = 45.0_lkind*dtr;  if (present(olatd)) olat=olatd*dtr
  wlon = 15.0_lkind*dtr;  if (present(wlond)) wlon=wlond*dtr
  wlat = 15.0_lkind*dtr;  if (present(wlatd)) wlat=wlatd*dtr
  rlon = 0.0_lkind     ;  if (present(rlond)) rlon=rlond*dtr
  rlat = 0.0_lkind     ;  if (present(rlatd)) rlat=rlatd*dtr

! compute gaussian-shaped mountain
    do j=1,size(lat(:))
      dy = abs(lat(j) - real(olat, FMS_TOP_KIND_))   ! dist from y origin
      yy = max(0.0_lkind, dy-real(rlat, FMS_TOP_KIND_))/real(wlat, FMS_TOP_KIND_)
      do i=1,size(lon(:))
        dx = abs(lon(i) - real(olon, FMS_TOP_KIND_)) ! dist from x origin
        dx = min(dx, abs(dx-tpi))  ! To ensure that: -pi <= dx <= pi
        xx = max(0.0_lkind, dx-real(rlon, FMS_TOP_KIND_))/real(wlon, FMS_TOP_KIND_)
        zsurf(i,j) = real(height, FMS_TOP_KIND_)*exp(-xx**2 - yy**2)
      enddo
    enddo

end function GET_GAUSSIAN_TOPOG_

!> @}
! close documentation grouping
