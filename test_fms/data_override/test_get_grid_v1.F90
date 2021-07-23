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

!> @brief  This programs tests calls to get_grid_version_1 used by data_override

program test_get_grid_v1

use netcdf,          only: nf90_create, nf90_clobber, nf90_64bit_offset, nf90_double, &
                           nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, &
                           nf90_close
use mpp_mod,         only: mpp_init, mpp_exit, mpp_root_pe, mpp_pe, mpp_error, FATAL, mpp_sync
use mpp_domains_mod, only: mpp_define_domains, mpp_define_io_domain, mpp_get_compute_domain, &
                           domain2d
use fms2_io_mod,     only: fms2_io_init
use platform_mod

use get_grid_version_mod, only : get_grid_version_1, deg_to_radian

implicit none

type(domain2d)                                   :: Domain !< 2D domain
integer                                          :: is, ie, js, je !< Starting and ending compute
                                                                   !! domain indices
integer                                          :: nlon, nlat !< Number of lat, lon in grid
real                                             :: min_lon, max_lon !< Maximum lat and lon
real, dimension(:,:), allocatable                :: lon, lat !< Lat and lon
integer                                          :: ncid, err !< Netcdf integers
integer                                          :: dimid1, dimid2, dimid3, dimid4 !< Dimensions IDs
integer                                          :: varid1, varid2, varid3, varid4, varid5 !< Variable IDs
real                                             :: lat_in(1), lon_in(1) !< Lat and lon to be written to file
real, dimension(:,:,:), allocatable              :: lat_vert_in, lon_vert_in !<Lat and lon vertices


lat_in = real(55.5, kind=r8_kind)
lon_in = real(44.5, kind=r8_kind)

call mpp_init
call fms2_io_init

!< Create a grid_spec file with information you need
if (mpp_pe() .eq. mpp_root_pe()) then
   err = nf90_create("grid_spec.nc", ior(nf90_clobber, nf90_64bit_offset), ncid)
   err = nf90_def_dim(ncid, 'xta', 1, dimid1)
   err = nf90_def_var(ncid, 'xta', nf90_double, (/dimid1/), varid1)

   err = nf90_def_dim(ncid, 'yta', 1, dimid2)
   err = nf90_def_var(ncid, 'yta', nf90_double, (/dimid2/), varid2)

   err = nf90_def_dim(ncid, 'gridlon_vert_t', 2, dimid3)
   err = nf90_def_dim(ncid, 'gridlat_vert_t', 2, dimid4)

   err = nf90_def_var(ncid, 'geolon_vert_t', nf90_double, (/dimid3,dimid4/), varid3)
   err = nf90_def_var(ncid, 'geolat_vert_t', nf90_double, (/dimid3,dimid4/), varid4)
   err = nf90_def_var(ncid, 'geolon_t',  nf90_double, (/dimid3,dimid4/), varid5)

   err = nf90_enddef(ncid)
   err = nf90_put_var(ncid, varid1, lon_in(1))
   err = nf90_put_var(ncid, varid2, lat_in(1))
   err = nf90_close(ncid)
endif

call mpp_sync()

nlon = 1
nlat = 1

!< Create a domain
call mpp_define_domains( (/1,nlon,1,nlat/), (/1, 1/), Domain, name='Atm')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_compute_domain(Domain,is,ie,js,je)

!< Call "get_grid_version_1" on a "atm" grid
allocate(lon(is:ie,js:je), lat(is:ie,js:je))
call get_grid_version_1("grid_spec.nc", "atm", Domain, is, ie, js, je, lon, lat, &
                        min_lon, max_lon)

!< Error checking:
if (lon(1,1) .ne. lon_in(1)*deg_to_radian) call mpp_error(FATAL,'test_get_grid_v1: lon is not the expected result')
if (lat(1,1) .ne. lat_in(1)*deg_to_radian) call mpp_error(FATAL,'test_get_grid_v1: lat is not the expected result')

!< Try again with ocean
lat = 0.
lon = 0.

call get_grid_version_1("grid_spec.nc", "ocn", Domain, is, ie, js, je, lon, lat, &
                        min_lon, max_lon)

!< Try again with ocean, "new_grid"
allocate(lat_vert_in(1,1,4), lon_vert_in(1,1,4))
lat_vert_in = real(55.5, kind=r8_kind)
lon_vert_in = real(65.5, kind=r8_kind)

!< Create a new grid_spec file with "x_T" instead of "geolon_t"
if (mpp_pe() .eq. mpp_root_pe()) then
   err = nf90_create("grid_spec.nc", ior(nf90_clobber, nf90_64bit_offset), ncid)
   err = nf90_def_dim(ncid, 'nlat', 1, dimid1)
   err = nf90_def_dim(ncid, 'nlon', 1, dimid2)
   err = nf90_def_dim(ncid, 'nz', 4, dimid3)

   err = nf90_def_var(ncid, 'x_T', nf90_double, (/dimid1, dimid2, dimid3/), varid1)
   err = nf90_def_var(ncid, 'x_vert_T', nf90_double, (/dimid1, dimid2, dimid3/), varid2)
   err = nf90_def_var(ncid, 'y_vert_T', nf90_double, (/dimid1, dimid2, dimid3/), varid3)

   err = nf90_enddef(ncid)
   err = nf90_put_var(ncid, varid2, lon_vert_in)
   err = nf90_put_var(ncid, varid3, lat_vert_in)

   err = nf90_close(ncid)
endif
call mpp_sync()

call get_grid_version_1("grid_spec.nc", "ocn", Domain, is, ie, js, je, lon, lat, &
                        min_lon, max_lon)

!< Error checking:
if (lon(1,1) .ne. sum(lon_vert_in)/4*deg_to_radian ) then
     call mpp_error(FATAL,'test_get_grid_v1: ocn, new grid, lon is not the expected result')
endif

if (lat(1,1) .ne. sum(lat_vert_in)/4*deg_to_radian ) then
     call mpp_error(FATAL,'test_get_grid_v1: ocn, new grid, lat is not the expected result')
endif

deallocate(lat_vert_in, lon_vert_in, lat, lon)

call mpp_exit()

end program test_get_grid_v1
