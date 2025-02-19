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

program test_collective_io
  use   fms_mod,         only: fms_init, fms_end, check_nml_error
  use   mpp_domains_mod, only: mpp_domains_set_stack_size, mpp_define_domains, mpp_define_io_domain, &
                               mpp_get_compute_domain, mpp_get_data_domain, domain2d, EAST, NORTH, CENTER, &
                               mpp_get_domain_tile_commid
  use   mpp_mod,         only: mpp_chksum, mpp_pe, mpp_root_pe, mpp_error, FATAL, input_nml_file, mpp_sync
  use   fms2_io_mod,     only: open_file, register_axis, register_variable_attribute, close_file, &
                               FmsNetcdfDomainFile_t, write_data, register_field, read_data

  implicit none

  !< Namelist variables
  integer, dimension(2) :: layout = (/1,6/)           !< Layout of the domain
  integer, dimension(2) :: io_layout = (/1,1/)        !< Io layout (currently overwritten to be the same as the
                                                      !! layout for reading, and 1,1 for writting)
  integer               :: nx = 96                    !< Size of the "x" dimension
  integer               :: ny = 96                    !< Size of the "y" dimension
  integer               :: nz = 65                    !< Size of the "z" dimension

  !< Other local variables
  type(FmsNetcdfDomainFile_t)           :: fileobj      !< FMS2_io fileobj
  type(domain2d)                        :: Domain_write !< Domain of the data for when writing the file
                                                        !! (it uses an io_layout of 1,1)
  type(domain2d)                        :: Domain_read  !< Domain of the data for when reading the file
                                                        !! (it uses an io_layout equal to the layout)
  character(len=6), dimension(3)        :: names        !< Dimension names for the dummy variables
  real, allocatable, dimension(:,:,:)   :: sst_in       !< Data to be written
  real, allocatable, dimension(:,:,:)   :: sst_out      !< Data read in
  real, allocatable, dimension(:,:)     :: sst2d_in       !< Data to be written
  real, allocatable, dimension(:,:)     :: sst2d_out      !< Data read in
  integer                               :: io_status    !< Status after reading the namelist

  integer                               :: is, ie, js, je !< Starting and ending indices
  integer                               :: i, j, k

  namelist / test_collective_io_nml / layout, io_layout, nx, ny, nz

  call fms_init()
  read (input_nml_file, test_collective_io_nml, iostat=io_status)
  if (io_status > 0) call mpp_error(FATAL,'=>test_collective_io_nml: Error reading input.nml')

  call mpp_domains_set_stack_size(17280000)
  call mpp_define_domains( (/1,nx,1,ny/), layout, Domain_write)
  call mpp_define_io_domain(Domain_write, (/1,1/))
  call mpp_get_compute_domain(Domain_write, is, ie, js, je)

  io_layout = layout !TODO remove when io_layout that are not equal to the layout works
  call mpp_define_domains( (/1,nx,1,ny/), layout, Domain_read)
  call mpp_define_io_domain(Domain_read, io_layout)

  names(1) = "lon"
  names(2) = "lat"
  names(3) = "level"

  allocate(sst_in(is:ie, js:je, nz))
  allocate(sst2d_in(is:ie, js:je))
  do i = is, ie
    do j = js, je
      do k = 1, nz
        sst_in(i,j,k) = i*10000. + j + k/100.
      enddo
    enddo
  enddo
  sst2d_in(:,:) = sst_in(:,:,1)

  allocate(sst_out(is:ie, js:je, nz))
  allocate(sst2d_out(is:ie, js:je))
  sst_out = real(-999.99)
  sst2d_out = real(-999.99)

  if (open_file(fileobj, "test_collective_io.nc", "overwrite", Domain_write, nc_format="netcdf4")) then
    call register_axis(fileobj, names(1), "x")
    call register_axis(fileobj, names(2), "y")
    call register_axis(fileobj, names(3), nz)

    call register_field(fileobj, "sst_3d", "double", names(1:3))
    call write_data(fileobj, "sst_3d", sst_in)

    call register_field(fileobj, "sst_2d", "double", names(1:2))
    call write_data(fileobj, "sst_2d", sst2d_in)

    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Unable to open the file for reading")
  endif
  call mpp_sync()

  fileobj%use_collective = .true.
  fileobj%tile_comm = mpp_get_domain_tile_commid(Domain_read)
  if (open_file(fileobj, "test_collective_io.nc", "read", Domain_read, nc_format="netcdf4")) then
    names(1) = "lon"
    names(2) = "lat"
    call register_axis(fileobj, "lon", "x")
    call register_axis(fileobj, "lat", "y")

    call read_data(fileobj, "sst_3d", sst_out)
    call read_data(fileobj, "sst_2d", sst2d_out)
    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Unable to open the file for reading")
  endif
  call mpp_sync()

  if (mpp_chksum(sst_in) .ne. mpp_chksum(sst_out)) &
    call mpp_error(FATAL, "Checksums don't match for the 3d variable")

  if (mpp_chksum(sst_in) .ne. mpp_chksum(sst_out)) &
    call mpp_error(FATAL, "Checksums don't match for the 3d variable")

  call fms_end()
end program