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

program test_io_with_mask

!> @brief  This programs tests fms2io/include/domain_write ability to write
!! data when the domain contains a mask table. For the points that are
!! masked out, no data should be writen.

use   mpp_domains_mod
use   mpp_mod
use   data_override_mod
use   fms2_io_mod
use   time_manager_mod, only: set_calendar_type, time_type, set_date, NOLEAP
use,  intrinsic :: iso_fortran_env, only : real64
use   mpi,             only: mpi_barrier, mpi_comm_world
use   netcdf

implicit none

integer, dimension(2)                 :: layout = (/2,3/) !< Domain layout
integer                               :: nlon             !< Number of points in x axis
integer                               :: nlat             !< Number of points in y axis
type(domain2d)                        :: Domain           !< Domain with mask table
real(kind=real64), allocatable, dimension(:,:) :: runoff  !< Data to be written
integer                               :: is, isc          !< Starting x index
integer                               :: ie, iec          !< Ending x index
integer                               :: js, jsc          !< Starting y index
integer                               :: je, jec          !< Ending y index
type(time_type)                       :: Time             !< Time
integer                               :: i, j             !< Helper indices
integer                               :: ncid             !< Netcdf file id
integer                               :: err              !< Error Code
integer                               :: dim1d, dim2d, dim3d, dim4d    !< Dimension ids
integer                               :: varid, varid2, varid3, varid4 !< Variable ids
real(kind=real64), allocatable, dimension(:,:,:) :: runoff_in          !< Data to be written to file
real(kind=real64)                     :: expected_result  !< Expected result from data_override

call mpp_init
call fms2_io_init
call data_override_init

!< Create some files needed by data_override!
if (mpp_pe() .eq. mpp_root_pe()) then
   allocate(runoff_in(1440, 1080, 10))
   do i = 1, 10
       runoff_in(:,:,i) = real(i, real64)
   enddo

   err = nf90_create('INPUT/grid_spec.nc', ior(nf90_clobber, nf90_64bit_offset), ncid)
   err = nf90_def_dim(ncid, 'str', 255, dim1d)
   err = nf90_def_var(ncid, 'ocn_mosaic_file', nf90_char, (/dim1d/), varid)
   err = nf90_enddef(ncid)
   err = nf90_put_var(ncid, varid, "ocean_mosaic.nc")
   err = nf90_close(ncid)

   err = nf90_create('INPUT/ocean_mosaic.nc', ior(nf90_clobber, nf90_64bit_offset), ncid)
   err = nf90_def_dim(ncid, 'str', 255, dim1d)
   err = nf90_def_dim(ncid, 'ntiles', 1, dim2d)
   err = nf90_def_var(ncid, 'gridfiles', nf90_char, (/dim1d, dim2d/), varid)
   err = nf90_enddef(ncid)
   err = nf90_put_var(ncid, varid, "ocean_hgrid.nc")
   err = nf90_close(ncid)

   err = nf90_create('INPUT/ocean_hgrid.nc', ior(nf90_clobber, nf90_64bit_offset), ncid)
   err = nf90_def_dim(ncid, 'nx', 2880, dim1d)
   err = nf90_def_dim(ncid, 'ny', 2160, dim2d)
   err = nf90_def_dim(ncid, 'nxp', 2881, dim3d)
   err = nf90_def_dim(ncid, 'nyp', 2161, dim4d)
   err = nf90_def_var(ncid, 'x', nf90_double, (/dim3d, dim4d/), varid)
   err = nf90_def_var(ncid, 'y', nf90_double, (/dim3d, dim4d/), varid2)
   err = nf90_def_var(ncid, 'area', nf90_double, (/dim1d, dim2d/), varid3)
   err = nf90_enddef(ncid)
   err = nf90_close(ncid)

   err = nf90_create('INPUT/runoff.daitren.clim.1440x1080.v20180328.nc', ior(nf90_clobber, nf90_64bit_offset), ncid)
   err = nf90_def_dim(ncid, 'i', 1440, dim1d)
   err = nf90_def_dim(ncid, 'j', 1080, dim2d)
   err = nf90_def_dim(ncid, 'time', nf90_unlimited, dim3d)

   err = nf90_def_var(ncid, 'i', nf90_double, (/dim1d/), varid3)
   err = nf90_put_att(ncid, varid3, "cartesian_axis", "x")

   err = nf90_def_var(ncid, 'j', nf90_double, (/dim2d/), varid4)
   err = nf90_put_att(ncid, varid4, "cartesian_axis", "y")

   err = nf90_def_var(ncid, 'time', nf90_double, (/dim3d/), varid2)
   err = nf90_put_att(ncid, varid2,  "cartesian_axis", "T")
   err = nf90_put_att(ncid, varid2, "calendar", "noleap")
   err = nf90_put_att(ncid, varid2, "units", "days since 0001-01-01 00:00:00")
   err = nf90_def_var(ncid, 'runoff', nf90_double, (/dim1d, dim2d, dim3d/), varid)

   err = nf90_enddef(ncid)
   err = nf90_put_var(ncid, varid, runoff_in)
   err = nf90_put_var(ncid, varid2, (/1., 2., 3., 5., 6., 7., 8., 9., 10., 11./))
   err = nf90_close(ncid)

endif

!< Wait for the root PE to catch up
call mpi_barrier(mpi_comm_world, err)

!< This is the actual test code:

call set_calendar_type(NOLEAP)
Time = set_date(1,1,5,0,0,0)

nlon = 1440
nlat = 1080

!< Create a domain nlonXnlat with mask
call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, xhalo=2, yhalo=2, name='test_data_override_emc')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, is, ie, js, je)

!< Set up the data
allocate(runoff(is:ie,js:je))

runoff = 999.

!< Initiliaze data_override
call data_override_init(Ocean_domain_in=Domain)
call data_override('OCN','runoff',runoff, Time)

!< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
!! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
expected_result = (real(3, real64)+ real(4, real64))/2

!! Data is only expected to be overriden for the compute domain -not at the halos.
call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
do i = is, ie
   do j = js, je
      if (i < isc .or. i > iec .or. j < jsc .or. j > jec) then
         !< This is the result at the halos it should 999.
         if (runoff(i,j) .ne. 999.) then
            print *, "for i=", i, " and j=", j, " runoff=", runoff(i,j)
            call mpp_error(FATAL, "test_data_override_ongrid: Result is different from expected answer!")
         endif
      else
         if (runoff(i,j) .ne. expected_result) then
            print *, "for i=", i, " and j=", j, " runoff=", runoff(i,j)
            call mpp_error(FATAL, "test_data_override_ongrid: Result is different from expected answer!")
         endif
      endif
   enddo
enddo

deallocate(runoff)

call mpp_exit

end program test_io_with_mask
