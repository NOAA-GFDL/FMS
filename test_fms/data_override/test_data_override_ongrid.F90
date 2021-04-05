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

program test_data_override_ongrid

!> @brief  This programs tests data_override ability to override data for an
!! on grid case

use   mpp_domains_mod,   only: mpp_define_domains, mpp_define_io_domain, mpp_get_data_domain, &
                               mpp_domains_set_stack_size, mpp_get_compute_domain, domain2d
use   mpp_mod,           only: mpp_init, mpp_exit, mpp_pe, mpp_root_pe, mpp_error, FATAL, &
                               input_nml_file
use   data_override_mod, only: data_override_init, data_override
use   fms2_io_mod,       only: fms2_io_init
use   time_manager_mod,  only: set_calendar_type, time_type, set_date, NOLEAP
use   mpi,               only: mpi_barrier, mpi_comm_world
use   netcdf,            only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, &
                               nf90_close, nf90_put_att, nf90_clobber, nf90_64bit_offset, nf90_char, &
                               nf90_double, nf90_unlimited

implicit none

integer, dimension(2)                 :: layout = (/2,3/) !< Domain layout
integer                               :: nlon             !< Number of points in x axis
integer                               :: nlat             !< Number of points in y axis
type(domain2d)                        :: Domain           !< Domain with mask table
real, allocatable, dimension(:,:)     :: runoff           !< Data to be written
integer                               :: is               !< Starting x index
integer                               :: ie               !< Ending x index
integer                               :: js               !< Starting y index
integer                               :: je               !< Ending y index
type(time_type)                       :: Time             !< Time
integer                               :: i                !< Helper indices
integer                               :: ncid             !< Netcdf file id
integer                               :: err              !< Error Code
integer                               :: dim1d, dim2d, dim3d, dim4d    !< Dimension ids
integer                               :: varid, varid2, varid3, varid4 !< Variable ids
real, allocatable, dimension(:,:,:)   :: runoff_in        !< Data to be written to file
real                                  :: expected_result  !< Expected result from data_override
integer                               :: nhalox=2, nhaloy=2
integer                               :: io_status

namelist / test_data_override_ongrid_nml / nhalox, nhaloy

call mpp_init
call fms2_io_init

read (input_nml_file, test_data_override_ongrid_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>test_data_override_ongrid: Error reading input.nml')

!< Create some files needed by data_override!
if (mpp_pe() .eq. mpp_root_pe()) then
   allocate(runoff_in(1440, 1080, 10))
   do i = 1, 10
       runoff_in(:,:,i) = real(i)
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

   deallocate(runoff_in)

endif

!< Wait for the root PE to catch up
call mpi_barrier(mpi_comm_world, err)

!< This is the actual test code:

call set_calendar_type(NOLEAP)

nlon = 1440
nlat = 1080

!< Create a domain nlonXnlat with mask
call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, xhalo=nhalox, yhalo=nhaloy, name='test_data_override_emc')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, is, ie, js, je)

print *, nhalox, nhaloy

!< Set up the data
allocate(runoff(is:ie,js:je))

runoff = 999.

!< Initiliaze data_override
call data_override_init(Ocean_domain_in=Domain)

!< Run it when time=3
Time = set_date(1,1,4,0,0,0)
call data_override('OCN','runoff',runoff, Time)
!< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
!! equal to the data at time=3, which is 3.
expected_result = real(3.)
call compare_data(Domain, runoff, expected_result)

!< Run it when time=4
runoff = 999.
Time = set_date(1,1,5,0,0,0)
call data_override('OCN','runoff',runoff, Time)
!< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
!! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
expected_result = (real(3.)+ real(4.))/2
call compare_data(Domain, runoff, expected_result)

deallocate(runoff)

call mpp_exit

contains

subroutine compare_data(Domain, actual_result, expected_result)
type(domain2d), intent(in)            :: Domain           !< Domain with mask table
real, intent(in)                      :: expected_result  !< Expected result from data_override
real, dimension(:,:), intent(in)      :: actual_result    !< Result from data_override

integer                               :: xsizec, ysizec   !< Size of the compute domain
integer                               :: xsized, ysized   !< Size of the data domain
integer                               :: nx, ny           !< Size of acual_result
integer                               :: nhalox, nhaloy   !< Size of the halos
integer                               :: i, j             !< Helper indices

!< Data is only expected to be overriden for the compute domain -not at the halos.
call mpp_get_compute_domain(Domain, xsize=xsizec, ysize=ysizec)
call mpp_get_data_domain(Domain, xsize=xsized, ysize=ysized)

!< Note that actual_result has indices at (1:nx,1:ny) not (is:ie,js:je)
nhalox= (xsized-xsizec)/2
nhaloy = (ysized-ysizec)/2
nx = size(actual_result, 1)
ny = size(actual_result, 2)

do i = 1, nx
   do j = 1, ny
      if (i <= nhalox .or. i > (nx-nhalox) .or. j <= nhaloy .or. j > (ny-nhaloy)) then
         !< This is the result at the halos it should 999.
         if (actual_result(i,j) .ne. 999.) then
            print *, "for i=", i, " and j=", j, " result=", actual_result(i,j)
            call mpp_error(FATAL, "test_data_override_ongrid: Data was overriden in the halos!!")
         endif
      else
         if (actual_result(i,j) .ne. expected_result) then
            print *, "for i=", i, " and j=", j, " result=", actual_result(i,j)
            call mpp_error(FATAL, "test_data_override_ongrid: Result is different from expected answer!")
         endif
      endif
   enddo
enddo

end subroutine

end program test_data_override_ongrid
