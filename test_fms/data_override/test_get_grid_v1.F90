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

use get_grid_version_mod, only : get_grid_version_1
use constants_mod, only: DEG_TO_RAD

implicit none

integer, parameter                         :: lkind = TEST_FMS_KIND_ !< r4_kind or r8_kind
type(domain2d)                             :: Domain !< 2D domain
integer                                    :: nlon, nlat !< Number of lat, lon in grid
real(lkind)                                :: min_lon, max_lon !< Maximum lat and lon
real(lkind), dimension(:,:), allocatable   :: lon, lat !< Lat and lon
integer                                    :: ncid, err !< Netcdf integers
integer                                    :: dimid1, dimid2, dimid3, dimid4 !< Dimensions IDs
integer                                    :: varid1, varid2, varid3, varid4, varid5 !< Variable IDs
real(lkind)                                :: lat_in(1), lon_in(1) !< Lat and lon to be written to file
real(lkind), dimension(:,:,:), allocatable :: lat_vert_in, lon_vert_in !<Lat and lon vertices

lat_in = 55.5_lkind
lon_in = 44.5_lkind

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

!< Call "get_grid_version_1" on a "atm" grid
call get_grid_version_1("grid_spec.nc", "atm", Domain, lon, lat, min_lon, max_lon)

!< Error checking:
if (lon(1,1) .ne. lon_in(1)*real(DEG_TO_RAD, lkind)) &
  & call mpp_error(FATAL,'test_get_grid_v1: lon is not the expected result')
if (lat(1,1) .ne. lat_in(1)*real(DEG_TO_RAD, lkind)) &
  & call mpp_error(FATAL,'test_get_grid_v1: lat is not the expected result')

call get_grid_version_1("grid_spec.nc", "ocn", Domain, lon, lat, min_lon, max_lon)

!< Try again with ocean, "new_grid"
allocate(lat_vert_in(1,1,4), lon_vert_in(1,1,4))
lat_vert_in = 55.5_lkind
lon_vert_in = 65.5_lkind

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

call get_grid_version_1("grid_spec.nc", "ocn", Domain, lon, lat, min_lon, max_lon)

!< Error checking:
if (lon(1,1) .ne. sum(lon_vert_in)/4._lkind * real(DEG_TO_RAD, lkind) ) then
     call mpp_error(FATAL,'test_get_grid_v1: ocn, new grid, lon is not the expected result')
endif

if (lat(1,1) .ne. sum(lat_vert_in)/4._lkind * real(DEG_TO_RAD, lkind) ) then
     call mpp_error(FATAL,'test_get_grid_v1: ocn, new grid, lat is not the expected result')
endif

deallocate(lat_vert_in, lon_vert_in, lat, lon)

call mpp_exit()

end program test_get_grid_v1
