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

subroutine atmosphere_restart_file(domain, nz, nt, debug_flag)

  type(domain2d), intent(in), target :: domain
  integer, intent(in) :: nz
  integer, intent(in) :: nt
  logical, intent(in), optional :: debug_flag

  type(FmsNetcdfDomainFile_t) :: fileobj
  type(domain2d), pointer :: io_domain
  integer :: isc
  integer :: nx
  integer :: jsc
  integer :: ny
  integer :: isd
  integer :: nxd
  integer :: jsd
  integer :: nyd
  integer :: ioff
  integer :: joff
  real(kind=r8_kind), dimension(:), allocatable :: lon_data
  real(kind=r8_kind), dimension(:), allocatable :: lat_data
  real(kind=r8_kind), dimension(:), allocatable :: lev_data
  real(kind=r8_kind), dimension(:,:,:), allocatable :: temperature_data
  real(kind=r8_kind), dimension(:,:,:), allocatable :: temperature_data_with_halos
  real(kind=r8_kind), dimension(:), allocatable :: avg_t_data
  integer(kind=i4_kind), dimension(:,:), allocatable :: surf_t_data
  integer :: ndims
  integer, dimension(:), allocatable :: dim_sizes
  integer(kind=i8_kind), dimension(3) :: out_chksum
  integer(kind=i8_kind), dimension(3) :: in_chksum
  integer, dimension(3) :: chksum_unlim_level
  character(len=256) :: att
  character(len=256), dimension(:), allocatable :: buffer
  character(len=256), dimension(:), allocatable :: buffer2
  character(len=256) :: str_buffer
  character(len=256) :: str_buffer2
  logical :: debug
  integer :: err
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  character(len=32) :: tstamp

  !Initialize.
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  if (present(debug_flag)) then
    debug = debug_flag
  else
    debug = .false.
  endif
  call random_seed()

  if (debug) then
    !Tell which test is running.
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
  call mpp_get_compute_domain(io_domain, xbegin=isc, xsize=nx, ybegin=jsc, ysize=ny)
  call mpp_get_data_domain(io_domain, xbegin=isd, xsize=nxd, ybegin=jsd, ysize=nyd)

  !Write a restart file.
  !Open a restart file and initialize the file object.
  call open_check(open_file(fileobj, "atmosphere.nc", "overwrite", &
                            domain, nc_format="64bit", is_restart=.true.))

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

  !Longitude2 (domain "x" dimension).
  call register_axis(fileobj, "lon2", "x")
  call register_field(fileobj, "lon2", "double", (/"lon2"/))
  do i = 1, nx
    lon_data(i) = 2.*3.14*lon_data(i)/(360.)
  enddo
  call write_data(fileobj, "lon2", lon_data)
  call register_variable_attribute(fileobj, "lon2", "units", "radians", str_len=len_trim("radians"))

  !Latitude (domain "y" dimension).
  call register_axis(fileobj, "lat", "y")
  call register_field(fileobj, "lat", "double", (/"lat"/))
  allocate(lat_data(ny))
  do i = 1, ny
    lat_data(i) = 65. + 0.05*i
  enddo
  call write_data(fileobj, "lat", lat_data)
  call register_variable_attribute(fileobj, "lat", "units", "degrees", str_len=len_trim("degrees"))

  !Latitude2 (domain "y" dimension).
  call register_axis(fileobj, "lat2", "y")
  call register_field(fileobj, "lat2", "double", (/"lat2"/))
  do i = 1, ny
    lat_data(i) = 2.*3.14*(lat_data(i))/360.
  enddo
  call write_data(fileobj, "lat2", lat_data)
  call register_variable_attribute(fileobj, "lat", "units2", "radians", str_len=len_trim("radians"))

  !Height.
  call register_axis(fileobj, "lev", nz)
  call register_field(fileobj, "lev", "double", (/"lev"/))
  allocate(lev_data(nz))
  do i = 1, nz
    lev_data(i) = 0.1 + 20.2*i
  enddo
  call write_data(fileobj, "lev", lev_data)
  call register_variable_attribute(fileobj, "lev", "units", "mb", str_len=len_trim("mb"))

  !Height2.
  call register_axis(fileobj, "lay", nz-1)
  call register_field(fileobj, "lay", "double", (/"lay"/))
  do i = 1, nz-1
    lev_data(i) = lev_data(i) - 0.5
  enddo
  call write_data(fileobj, "lay", lev_data(1:nz-1))
  call register_variable_attribute(fileobj, "lay", "units", "mb", str_len=len_trim("mb"))

  !Time.
  call register_axis(fileobj, "time", unlimited)
  call register_field(fileobj, "time", "float", (/"time"/))
  call register_variable_attribute(fileobj, "time", "units", "years", str_len=len_trim("years"))

  !String length.
  call register_axis(fileobj, "strlen", 256)

  !Add a restart variable.
  allocate(temperature_data(nx, ny, nz))
  temperature_data(:,:,:) = 0.0
  call register_restart_field(fileobj, "temperature", temperature_data, &
                              (/"lon ", "lat ", "lev ", "time"/))
  call register_variable_attribute(fileobj, "temperature", "units", "K", str_len=len_trim("k"))

  !Add a non-domain decomposed variable.
  allocate(avg_t_data(nz-1))
  call register_field(fileobj, "avg_layer_temperature", "double", (/"lay ", "time"/))
  call register_variable_attribute(fileobj, "avg_layer_temperature", "units", "K", str_len=len_trim("k"))

  !Add a domain decomposed variable.
  allocate(surf_t_data(nx, ny))
  call register_field(fileobj, "surface_temperature", "int", (/"lon2", "lat2", "time"/))
  call register_variable_attribute(fileobj, "surface_temperature", "units", "K", str_len=len_trim("k"))

  !Add scalar and 1D string variables.
  call register_field(fileobj, "filename", "char", (/"strlen"/))
  call register_field(fileobj, "strings", "char", (/"strlen", "lev   "/))
  call write_data(fileobj, "filename", "file1.nc")
  allocate(buffer(nz))
  do i = 1, nz
    buffer(i) = ""
    write(buffer(i), *) i, "string"
    buffer(i) = trim(adjustl(buffer(i)))
  enddo
  call write_data(fileobj, "strings", buffer)

  !Add a variable whose user-level buffers include space for halos.
  allocate(temperature_data_with_halos(nxd, nyd, nz))
  temperature_data(:,:,:) = 0.0
  call register_restart_field(fileobj, "temperature_ignore_halos", &
                              temperature_data_with_halos, &
                              (/"lon ", "lat ", "lev ", "time"/))

  !Initialize some random temperature data.
  do j = 1, nz
    do k = 1, ny
      do l = 1, nx
        call random_number(temperature_data(l, k, j))
        temperature_data(l, k, j) = temperature_data(l, k, j) + 300.05
      enddo
    enddo
  enddo
  do i = 1, nz-1
    avg_t_data(i) = sum(temperature_data(:, :, i))/(nx*ny)
  enddo
  do i = 1, ny
    do j = 1, nx
      surf_t_data(j, i) = temperature_data(j, i, 1)
    enddo
  enddo

  !Perform a "simulation".
  do i = 1, nt
    do j = 1, nz
      do k = 1, ny
        do l = 1, nx
          temperature_data(l, k, j) = temperature_data(l, k, j) + 0.5*mpp_pe()
        enddo
      enddo
    enddo
    do j = 1, nz-1
      avg_t_data(j) = sum(temperature_data(:, :, j))/(nx*ny)
    enddo
    call write_data(fileobj, "avg_layer_temperature", avg_t_data, unlim_dim_level=i)
    if (i .eq. 1) then
      out_chksum(2) = mpp_chksum(avg_t_data, pelist=(/mpp_pe()/))
      chksum_unlim_level(2) = i
    endif
    do j = 1, ny
      do k = 1, nx
        surf_t_data(k, j) = int(temperature_data(k, j, 1))
      enddo
    enddo
    call write_data(fileobj, "surface_temperature", surf_t_data, unlim_dim_level=i)
    if (i .eq. 2) then
      out_chksum(3) = mpp_chksum(surf_t_data, pelist=(/mpp_pe()/))
      chksum_unlim_level(3) = i
    endif
    ioff = isc-isd
    joff = jsc-jsd
    do j = 1, nz
      do k = 1, ny
        do l = 1, nx
          temperature_data_with_halos(l+ioff,k+joff,j) = temperature_data(l,k,j)
        enddo
      enddo
    enddo

    !Write out the times of each "time step".
    call write_data(fileobj, "time", 1950. + real(i,kind=r4_kind), unlim_dim_level=i)
  enddo

  !Write out the restart file.
  call write_restart(fileobj, unlim_dim_level=nt)
  out_chksum(1) = mpp_chksum(temperature_data, pelist=(/mpp_pe()/))
  chksum_unlim_level(1) = nt

  !Close the file and destroy the file object.
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Prepare to read the data back in.
  deallocate(lon_data)
  deallocate(lat_data)
  deallocate(lev_data)
  deallocate(temperature_data)
  deallocate(avg_t_data)
  deallocate(surf_t_data)
  deallocate(temperature_data_with_halos)

  !Check if a non-existent file exists (just to see if this feature
  !works.
  if (open_file(fileobj, "atmosphere.foobar.nc", "read", domain, &
                nc_format="64bit", is_restart=.true.)) then
    call mpp_error(FATAL, "Found non-existent file.")
  endif

  !Re-open the restart file and re-initialize the file object.
  call open_check(open_file(fileobj, "atmosphere.nc", "read", domain, &
                            nc_format="64bit", is_restart=.true.))

  !Associate the lon and lat dimensions with the
  !"x" and "y" dimensions of the domain.
  call register_axis(fileobj, "lon", "x")
  call register_axis(fileobj, "lat", "y")
  call register_axis(fileobj, "lon2", "x")
  call register_axis(fileobj, "lat2", "y")

  !Read in and check the string variables.
  call get_dimension_size(fileobj, "lev", ndims)
  allocate(buffer2(ndims))
  call read_data(fileobj, "strings", buffer2)
  do i = 1, ndims
    if (trim(buffer2(i)) .ne. trim(buffer(i))) then
      call mpp_error(fatal, "Did not read in strings correctly.")
    endif
  enddo
  call read_data(fileobj, "strings", str_buffer, corner=ndims)
  if (trim(str_buffer) .ne. trim(buffer(ndims))) then
    call mpp_error(fatal, "Scalar from array of strings read failed.")
  endif
  call read_data(fileobj, "filename", str_buffer2)
  if (trim(str_buffer2) .ne. "file1.nc") then
    call mpp_error(fatal, "Did not read in filename variable correctly.")
  endif
  deallocate(buffer)
  deallocate(buffer2)

  !Allocate buffers associated with the temperature variables.
  ndims = get_variable_num_dimensions(fileobj, "temperature", broadcast=.true.)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, "temperature", dim_sizes, broadcast=.true.)
  allocate(temperature_data(nx, ny, dim_sizes(3)))
  temperature_data(:,:,:) = 0.0
  allocate(temperature_data_with_halos(nxd, nyd, dim_sizes(3)))
  temperature_data_with_halos(:,:,:) = 0.0
  deallocate(dim_sizes)

  ndims = get_variable_num_dimensions(fileobj, "avg_layer_temperature", broadcast=.true.)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, "avg_layer_temperature", dim_sizes, broadcast=.true.)
  allocate(avg_t_data(dim_sizes(1)))
  deallocate(dim_sizes)
  avg_t_data(:) = 0.
  allocate(surf_t_data(nx, ny))
  surf_t_data(:, :) = 0

  !Register the restart variable (so a pointer can point at the newly
  !allocated buffer.
  call register_restart_field(fileobj, "temperature", temperature_data)
  call register_restart_field(fileobj, "temperature_ignore_halos", &
                              temperature_data_with_halos)

  !Read in variable attribute.
  call get_variable_attribute(fileobj, "temperature", "units", att)
  if (trim(att) .ne. "K") then
    call mpp_error(fatal, "failed to read in string variable attribute.")
  endif

  !Read in the restart variable and checksum.
  call read_restart(fileobj, unlim_dim_level=chksum_unlim_level(1))
  in_chksum(1) = mpp_chksum(temperature_data, pelist=(/mpp_pe()/))
  call chksum_match(out_chksum(1), in_chksum(1), "temperature", debug)

  !Read in the non-restart variables and checksum.
  call chksum_match(out_chksum(1), in_chksum(1), "temperature_ignore_halos", debug)
  call read_data(fileobj, "avg_layer_temperature", avg_t_data, &
                 unlim_dim_level=chksum_unlim_level(2))
! in_chksum(2) = mpp_chksum(avg_t_data, pelist=(/mpp_pe()/))
! call chksum_match(out_chksum(2), in_chksum(2), "avg_layer_temperature", debug)
  call read_data(fileobj, "surface_temperature", surf_t_data, &
                 unlim_dim_level=chksum_unlim_level(3))
  in_chksum(3) = mpp_chksum(surf_t_data, pelist=(/mpp_pe()/))
  call chksum_match(out_chksum(3), in_chksum(3), "surface_temperature", debug)

  !Close the file.
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  if (mpp_pe() .eq. 0) then
    write(output_unit,'(a/)') &
      "Atmosphere (6-tile domain decomposed) file test ... " &
      //trim(green)//"passed."//trim(color_end)
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
end subroutine atmosphere_restart_file
