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

! This is a very simple test of FMS IO.

! Ed Hartnett, 6/10/20

program test_io_simple
  use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64, error_unit, output_unit
  use mpi
  use fms2_io_mod
  use netcdf_io_mod
  use mpp_domains_mod
  use mpp_mod
  use setup
  use netcdf
  
  type(Params) :: test_params
  type(domain2d) :: domain
  integer :: err
  type(domain2d), pointer :: io_domain
  type(FmsNetcdfDomainFile_t) :: fileobj
  integer, parameter :: ntiles = 6
  integer :: ncchksz = 64*1024 
  character (len = 10) :: netcdf_default_format = "64bit"
  integer :: header_buffer_val = 16384
  character (len = 100) :: testfile
  integer :: ncid
  integer :: numfilesatt

  ! Initialize.
  call init(test_params, ntiles)
  call create_cubed_sphere_domain(test_params, domain, (/1, 1/))
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  if (test_params%debug) then
     if (mpp_pe() .eq. 0) then
        write(error_unit,'(/a)') &
             "Running atmosphere (6-tile domain decomposed) restart file test ... "
     endif
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Get the sizes of the I/O compute and data domains.
  io_domain => mpp_get_io_domain(domain)
  if (.not. associated(io_domain)) then
     call mpp_error(fatal, "I/O domain is not associated.")
  endif
  call mpp_get_compute_domain(io_domain, xbegin=isc, xend=iec, xsize=nx, ybegin=jsc, yend=jec, ysize=ny)
  call mpp_get_data_domain(io_domain, xbegin=isd, xsize=nxd, ybegin=jsd, ysize=nyd)
  call mpp_get_compute_domain(io_domain, xbegin=isc_east, xend=iec_east, xsize=nx_east, position=east)
  call mpp_get_compute_domain(io_domain, ybegin=jsc_north, yend=jec_north, ysize=ny_north, position=north)

  call netcdf_io_init(ncchksz, header_buffer_val, netcdf_default_format)

  ! Open a restart file and initialize the file object.
  call open_check(open_file(fileobj, "test_io_simple.nc", "overwrite", &
       domain, nc_format="64bit", is_restart=.false.))

  ! Close the file.
  call close_file(fileobj)

  ! Check for expected netcdf file.
  if (mpp_pe() .eq. 0) then
     do i = 1, 6
        write(testfile,'(a,i1,a)') 'test_io_simple.tile', i, '.nc'
        err = nf90_open(testfile, nf90_nowrite, ncid)
        if (err .ne. 0) stop 1
        err = nf90_get_att(ncid, NF_GLOBAL, 'NumFilesInSet', numfilesatt)
        if (err .ne. 0) stop 10
        if (numfilesatt .ne. 1) stop 11
        err = nf90_close(ncid)
        if (err .ne. 0) stop 90
     end do
  endif

  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  call cleanup(test_params)
end program test_io_simple
