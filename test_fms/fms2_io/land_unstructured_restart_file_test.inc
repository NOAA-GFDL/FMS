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

subroutine land_unstructured_restart_file(domain, nz, nt, debug_flag)

  type(domainug), intent(in), target :: domain
  integer, intent(in) :: nz
  integer, intent(in) :: nt
  logical, intent(in), optional :: debug_flag

  type(FmsNetcdfUnstructuredDomainFile_t) :: fileobj
  type(domainug), pointer :: io_domain
  integer :: total_nx
  integer :: s
  integer, dimension(:), allocatable :: s2
  integer :: nx
  real(kind=real64), dimension(:), allocatable :: tile_data
  real(kind=real64), dimension(:), allocatable :: lev_data
  real(kind=real64), dimension(:,:), allocatable :: temperature_data
  integer(kind=int64) :: out_chksum
  integer :: ndims
  integer, dimension(:), allocatable :: dim_sizes
  integer(kind=int64) :: in_chksum
  logical :: debug
  integer :: err
  integer :: i
  integer :: j
  integer :: k

  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  if (mpp_pe() .eq. 0) then
    write(output_unit,'(/a)',advance="no") &
      "Running land (unstructured domain decomposed) restart file test ... "
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  if (present(debug_flag)) then
    debug = debug_flag
  else
    debug = .false.
  endif
  call random_seed()

  !Get the sizes of the I/O compute domain, allocate buffer,
  !and initialize some data.
  io_domain => mpp_get_ug_io_domain(domain)
  if (.not. associated(io_domain)) then
    call mpp_error(fatal, "I/O domain is not associated.")
  endif
  allocate(dim_sizes(mpp_get_UG_domain_npes(io_domain)))
  allocate(s2(mpp_get_UG_domain_npes(io_domain)))
  dim_sizes(:) = 0
  call mpp_get_ug_compute_domains(io_domain, begin=s2, size=dim_sizes)
  total_nx = sum(dim_sizes)
  deallocate(dim_sizes)
  call mpp_get_ug_compute_domain(io_domain, begin=s, size=nx)
  if (s2(1) .ne. 1) then
    s = s - s2(1) + 1
  endif
  deallocate(s2)
  allocate(temperature_data(nx, nz))
  do j = 1, nz
    do k = 1, nx
      call random_number(temperature_data(k, j))
      temperature_data(k, j) = temperature_data(k, j) + 300.05
    enddo
  enddo

  !Write a restart file.
  !Open a restart file and initialize the file object.
  call open_check(open_file(fileobj, "land_unstructured.nc", "overwrite", &
                            domain, nc_format="classic", is_restart=.true.))

  !Add some global attributes to the restart file.
  call register_global_attribute(fileobj, "globalatt1", real(7., kind=real64))
  call register_global_attribute(fileobj, "globalatt2", real(4., kind=real32))
  call register_global_attribute(fileobj, "globalatt3", int(3, kind=int32))
  call register_global_attribute(fileobj, "globalatt4", int(2, kind=int64))
  call register_global_attribute(fileobj, "globalatt5", "some text", str_len=len_trim("some text"))

  !Add dimensions and corresponding variables to the file.
  !Tile (unstructured domain "x" dimension).
  call register_axis(fileobj, "tile")
  call register_field(fileobj, "tile", "double", (/"tile"/))
  allocate(tile_data(nx))
  do i = 1, nx
    tile_data(i) = 2.234 + 1.7567*i
  enddo
  call write_data(fileobj, "tile", tile_data)
  call register_variable_attribute(fileobj, "tile", "units", "tile value", str_len=len_trim("tile value"))

  !Height.
  call register_axis(fileobj, "lev", nz)
  call register_field(fileobj, "lev", "double", (/"lev"/))
  allocate(lev_data(nz))
  do i = 1, nz
    lev_data(i) = 0.1111 + 20.236*i
  enddo
  call write_data(fileobj, "lev", lev_data)
  call register_variable_attribute(fileobj, "lev", "units", "mb", str_len=len_trim("mb"))

  !Time.
  call register_axis(fileobj, "time", unlimited)
  call register_field(fileobj, "time", "float", (/"time"/))
  call register_variable_attribute(fileobj, "time", "units", "seconds", str_len=len_trim("seconds"))

  !Add a restart variable.
  call register_restart_field(fileobj, "temperature", temperature_data, &
                              (/"tile","lev ","time"/))
  call register_variable_attribute(fileobj, "temperature", "units", "K", str_len=len_trim("k"))

  !Initialize some random temperature data.

  !Perform a "simulation".
  do i = 1, nt
    do j = 1, nz
      do k = 1, nx
        temperature_data(k, j) = temperature_data(k, j) + 0.5*mpp_pe()
      enddo
    enddo
  enddo

  !Calculate a checksum, that will be used to check that the data is
  !read back in correctly.
  out_chksum = mpp_chksum(temperature_data, pelist=(/mpp_pe()/))

  !Write out the restart file.
  call write_restart(fileobj, unlim_dim_level=1)

  !Close the file and destroy the file object.
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Prepare to read the data back in.
  deallocate(tile_data)
  deallocate(lev_data)
  deallocate(temperature_data)

  !Read restart with an explicit axis register call.
  !Re-open the restart file and re-initialize the file object.
  call open_check(open_file(fileobj, "land_unstructured.nc", "read", domain, &
                            nc_format="classic", is_restart=.true.))

  !Associate the "tile" dimension with the unstructured
  !dimension of the unstructured domain.  This
  !is the "explicit" axis register call.
  call register_axis(fileobj, "tile")
  ndims = get_variable_num_dimensions(fileobj, "tile", broadcast=.true.)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, "tile", dim_sizes, .true.)
  if (dim_sizes(1) .ne. total_nx) then
    call mpp_error(fatal, "unstructured axis is not the correct size.")
  endif
  deallocate(dim_sizes)

  !Allocate a buffer and associate it with the
  !temperature variable.
  ndims = get_variable_num_dimensions(fileobj, "temperature", broadcast=.true.)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, "temperature", dim_sizes, broadcast=.true.)
  allocate(temperature_data(dim_sizes(1), dim_sizes(2)))
  deallocate(dim_sizes)
  temperature_data(:,:) = 0.

  !Read in the restart variable and close the file.
  call read_data(fileobj, "temperature", temperature_data, unlim_dim_level=1)
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Check to see if the read-in data matches the data that was
  !written out.
  in_chksum = mpp_chksum(temperature_data(s:s+nx-1,:), pelist=(/mpp_pe()/))
  call chksum_match(out_chksum, in_chksum, "temperature", debug)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  if (mpp_pe() .eq. 0) then
    write(output_unit,'(a/)') trim(green)//"passed."//trim(color_end)
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
end subroutine land_unstructured_restart_file
