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

subroutine ocean_restart_file(domain, nz, nt, debug_flag)

  type(domain2d), intent(in), target :: domain
  integer, intent(in) :: nz
  integer, intent(in) :: nt
  logical, intent(in), optional :: debug_flag

  type(FmsNetcdfDomainFile_t) :: fileobj
  type(domain2d), pointer :: io_domain
  integer :: nx
  integer :: ny
  real(kind=r8_kind), dimension(:), allocatable :: lon_data
  real(kind=r8_kind), dimension(:), allocatable :: lat_data
  real(kind=r8_kind), dimension(:), allocatable :: lev_data
  real(kind=r8_kind), dimension(:,:,:), allocatable :: temperature_data
  integer(kind=i8_kind) :: out_chksum
  integer :: ndims
  integer, dimension(:), allocatable :: dim_sizes
  integer(kind=i8_kind) :: in_chksum
  logical :: debug
  integer :: err
  integer :: i
  integer :: j
  integer :: k
  integer :: l

  !Tell which test is running.
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  if (mpp_pe() .eq. 0) then
    write(output_unit,'(/a)',advance="no") &
        "Running ocean (1-tile domain decomposed) restart file test ... "
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Initialize.
  if (present(debug_flag)) then
    debug = debug_flag
  else
    debug = .false.
  endif
  call random_seed()

  !Get the sizes of the I/O compute domain and allocate buffer.
  io_domain => mpp_get_io_domain(domain)
  if (.not. associated(io_domain)) then
    call mpp_error(fatal, "I/O domain is not associated.")
  endif
  call mpp_get_compute_domain(io_domain, xsize=nx, ysize=ny)
  allocate(temperature_data(nz, nx, ny))
  temperature_data(:,:,:) = 0.0

  !Write a restart file.
  !Open a restart file and initialize the file object.
  call open_check(open_file(fileobj, "ocean.nc", "overwrite", domain, &
                            nc_format="64bit", is_restart=.true.))

  !Add dimensions and corresponding variables to the file.
  !Longitude (domain "x" dimension).
  call register_axis(fileobj, "lon", "x")
  call register_field(fileobj, "lon", "double", (/"lon"/))
  allocate(lon_data(nx))
  do i = 1, nx
    lon_data(i) = 33.33 + 1.2*i
  enddo
  call write_data(fileobj, "lon", lon_data)
  call register_variable_attribute(fileobj, "lon", "units", "degrees", str_len=len_trim("degrees"))

  !Latitude (domain "y" dimension).
  call register_axis(fileobj, "lat", "y")
  call register_field(fileobj, "lat", "double", (/"lat"/))
  allocate(lat_data(ny))
  do i = 1, ny
    lat_data(i) = 65. + 0.05*i
  enddo
  call write_data(fileobj, "lat", lat_data)
  call register_variable_attribute(fileobj, "lat", "units", "degrees", str_len=len_trim("degrees"))

  !Height.
  call register_axis(fileobj, "lev", nz)
  call register_field(fileobj, "lev", "double", (/"lev"/))
  allocate(lev_data(nz))
  do i = 1, nz
    lev_data(i) = 0.1 + 20.2*i
  enddo
  call write_data(fileobj, "lev", lev_data)
  call register_variable_attribute(fileobj, "lev", "units", "mb", str_len=len_trim("mb"))

  !Add a restart variable.
  call register_restart_field(fileobj, "temperature", temperature_data, &
                              (/"lev ","lon ","lat "/))
  call register_variable_attribute(fileobj, "temperature", "units", "K", str_len=len_trim("K"))

  !Initialize some random temperature data.
  do j = 1, ny
    do k = 1, nx
      do l = 1, nz
        call random_number(temperature_data(l, k, j))
        temperature_data(l, k, j) = temperature_data(l, k, j) + 300.05
      enddo
    enddo
  enddo

  !Perform a "simulation".
  do i = 1, nt
    do j = 1, ny
      do k = 1, nx
        do l = 1, nz
          temperature_data(l, k, j) = temperature_data(l, k, j) + 0.5*mpp_pe()
        enddo
      enddo
    enddo
  enddo

  !Calculate a checksum, that will be used to check that the data is
  !read back in correctly.
  out_chksum = mpp_chksum(temperature_data, pelist=(/mpp_pe()/))

  !Write out the restart file.
  call write_restart(fileobj)

  !Close the file and destroy the file object.
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Prepare to read the data back in.
  deallocate(lon_data)
  deallocate(lat_data)
  deallocate(lev_data)
  deallocate(temperature_data)

  !Read restart with explicit axis register calls.
  !Re-open the restart file and re-initialize the file object.
  call open_check(open_file(fileobj, "ocean.nc", "read", domain, nc_format="64bit", &
                            is_restart=.true.))

  !Associate the lon and lat dimensions with the
  !"x" and "y" dimensions of the domain.  These
  !are the "explicit" axis register calls.
  call register_axis(fileobj, "lon", "x")
  call register_axis(fileobj, "lat", "y")

  !Allocate a buffer and associate it with the
  !temperature variable.
  ndims = get_variable_num_dimensions(fileobj, "temperature", broadcast=.true.)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, "temperature", dim_sizes, broadcast=.true.)
  allocate(temperature_data(dim_sizes(1), nx, ny))
  call register_restart_field(fileobj, "temperature", temperature_data)
  deallocate(dim_sizes)
  temperature_data(:, :, :) = 0.

  !Read in the restart variable and close the file.
  call read_restart(fileobj)
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Check to see if the read-in data matches the data that was
  !written out.
  in_chksum = mpp_chksum(temperature_data, pelist=(/mpp_pe()/))
  call chksum_match(out_chksum, in_chksum, "temperature", debug)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  if (mpp_pe() .eq. 0) then
    write(output_unit, '(a/)') trim(green)//"passed."//trim(color_end)
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
end subroutine ocean_restart_file
