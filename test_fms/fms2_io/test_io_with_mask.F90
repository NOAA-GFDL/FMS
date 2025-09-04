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

program test_io_with_mask
!> @brief  This programs tests fms2io/include/domain_write ability to write
!! data when the domain contains a mask table. For the points that are
!! masked out, no data should be writen.
!! It also tests fms2io/include/domain_read ability to read the data when
!! the domain contains a mask table. For this case the masked data should
!! not be read.

use   mpp_domains_mod, only: mpp_domains_set_stack_size, mpp_define_domains, mpp_define_io_domain, &
                             mpp_get_compute_domain,domain2d
use   mpp_mod,         only: mpp_pe, mpp_root_pe, mpp_error, FATAL
use   fms2_io_mod,     only: open_file, register_axis, register_variable_attribute, close_file, &
                             FmsNetcdfDomainFile_t, write_data, register_field, read_data, &
                             parse_mask_table
use   fms_mod,         only: fms_init, fms_end
use   netcdf,          only: nf90_open, nf90_get_var, nf90_nowrite, NF90_NOERR, nf90_get_var, &
                             nf90_close
use   mpi,             only: mpi_barrier, mpi_comm_world
use   platform_mod

implicit none

integer, dimension(2)                 :: layout = (/2,3/) !< Domain layout
integer                               :: nlon             !< Number of points in x axis
integer                               :: nlat             !< Number of points in y axis
type(domain2d)                        :: Domain           !< Domain with mask table
real, dimension(:), allocatable       :: x                !< x axis data
real, dimension(:), allocatable       :: y                !< y axis data
real(kind=r8_kind), allocatable, dimension(:,:) :: sst     !< Data to be written
real(kind=r8_kind), allocatable, dimension(:,:) :: sst_in  !< Buffer where data will be read with netcdf
real(kind=r8_kind), allocatable, dimension(:,:) :: sst_in2 !< Buffer where data will be read with fms2io
logical, allocatable, dimension(:,:)  :: parsed_mask      !< Parsed masked
character(len=6), dimension(2)        :: names            !< Dimensions names
type(FmsNetcdfDomainFile_t)           :: fileobj          !< fms2io fileobj for domain decomposed
integer                               :: err              !< Return code.
integer                               :: ncid             !< File ID for checking file.
integer                               :: i, j             !< Helper integers
integer                               :: is               !< Starting x index
integer                               :: ie               !< Ending x index
integer                               :: js               !< Starting y index
integer                               :: je               !< Ending y index

call fms_init

nlon = 60
nlat = 60

!< Parse the mask table
allocate(parsed_mask(layout(1), layout(2)))
call parse_mask_table("the_mask", parsed_mask, 'test_io_with_mask')

!< Create a domain nlonXnlat with mask
call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_io_with_mask', maskmap=parsed_mask)
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_compute_domain(Domain, is, ie, js, je)

!< Set up the data
allocate(x(nlon), y(nlat))
allocate(sst(is:ie,js:je))
allocate(sst_in2(is:ie,js:je))

do i=1,nlon
  x(i) = i
enddo
do j=1,nlat
  y(j) = j
enddo

sst = real(7., kind=r8_kind)

!< Open a netCDF file and initialize the file object.
if (open_file(fileobj, "test_io_with_mask.nc", "overwrite", domain)) then
    !< Register the axis
    names(1) = "lon"
    names(2) = "lat"
    call register_axis(fileobj, "lon", "x")
    call register_axis(fileobj, "lat", "y")

    !< Register the variable and Write out the data
    call register_field(fileobj, "sst", "double", names(1:2))
    call register_variable_attribute(fileobj, "sst", "_FillValue", real(999., kind=r8_kind))
    call write_data(fileobj, "sst", sst)

    !< Close the file
    call close_file(fileobj)
else
   call mpp_error(FATAL, "test_io_with_mask: error opening the file for writing")
endif

!< Read the file back to check if it was done right!
if (mpp_pe() .eq. mpp_root_pe()) then
   err = nf90_open("test_io_with_mask.nc", nf90_nowrite, ncid)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_io_with_mask: error opening the file for reading")

   allocate(sst_in(nlon,nlat))
   err = nf90_get_var(ncid, 1, sst_in)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_io_with_mask: error reading from the file")

   !< x: 1-30 y: 1-20 are masked out, so sst_in for this values has to be equal to the fill value
   !! For the other points the data should be equal to real(7., kind=r8_kind)

   do i=1,nlon
       do j=1,nlat
          if (i > 30 .or. j > 20) then
             if (sst_in(i,j) .ne. real(7., kind=r8_kind)) then
                 print *, 'i=', i, ' j=', j, ' sst_in=', sst_in(i,j)
                 call mpp_error(FATAL, "test_io_with_mask: the unmasked data read in is not correct")
             endif
          else
             if (sst_in(i,j) .ne. real(999., kind=r8_kind)) then
                 print *, 'i=', i, ' j=', j, ' sst_in=', sst_in(i,j)
                 call mpp_error(FATAL, "test_io_with_mask: the masked data read in is not correct")
             endif
          endif
       enddo
   enddo

   err = nf90_close(ncid)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_io_with_mask: error closing the file")
endif

!< Wait for the root pe to catch up!
call mpi_barrier(mpi_comm_world, err)

!< Read the file back using fms2io
sst_in2 = 0.
if (open_file(fileobj, "test_io_with_mask.nc", "read", domain)) then
   names(1) = "lon"
   names(2) = "lat"
   call register_axis(fileobj, "lon", "x")
   call register_axis(fileobj, "lat", "y")

   !< Register the variable and read out the data
   call register_field(fileobj, "sst", "double", names(1:2))

   call read_data(fileobj, "sst", sst_in2)
   call close_file(fileobj)

   do i=is,ie
       do j=js,je
             if (sst_in2(i,j) .ne. real(7., kind=r8_kind)) then
                 print *, 'i=', i, ' j=', j, ' sst_in=', sst_in2(i,j)
                 call mpp_error(FATAL, "test_io_with_mask: the unmasked data read in is not correct")
             endif
       enddo
   enddo

endif


call fms_end
end program test_io_with_mask
