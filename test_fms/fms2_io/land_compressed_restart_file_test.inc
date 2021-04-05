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

subroutine land_compressed_restart_file(nz, nt, debug_flag)

  integer, intent(in) :: nz
  integer, intent(in) :: nt
  logical, intent(in), optional :: debug_flag

  type(FmsNetcdfFile_t) :: fileobj
  integer :: my_rank
  integer :: num_ranks
  integer, dimension(:), allocatable :: pelist
  integer, dimension(:), allocatable :: npes_corner
  integer, dimension(:), allocatable :: npes_nelems
  integer :: nx
  integer :: total_nx
  real(kind=real64), dimension(:,:), allocatable :: temperature_data
  real(kind=real64), dimension(:), allocatable :: cdim_data
  real(kind=real64), dimension(:), allocatable :: lev_data
  integer(kind=int64) :: out_chksum
  integer :: ndims
  integer, dimension(:), allocatable :: dim_sizes
  integer(kind=int64) :: in_chksum
  logical :: debug
  integer :: err
  integer :: i
  integer :: j
  integer :: k
  integer :: s

  !Tell which test is running.
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  if (mpp_pe() .eq. 0) then
    write(output_unit,'(/a)',advance="no") &
      "Running land (compressed) restart file test ... "
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

  !Initalize some data.  Each rank will have (my_rank+1)*5
  !data elements.
  my_rank = mpp_pe()
  num_ranks = mpp_npes()
  allocate(pelist(num_ranks))
  allocate(npes_corner(num_ranks))
  allocate(npes_nelems(num_ranks))
  nx = 1
  do i = 1, npes
    pelist(i) = i - 1
    npes_corner(i) = nx
    npes_nelems(i) = i*5
    nx = nx + npes_nelems(i)
  enddo
  nx = npes_nelems(my_rank+1)
  total_nx = sum(npes_nelems)
  allocate(temperature_data(nx, nz))
  do i = 1, nz
    do j = 1, nx
      call random_number(temperature_data(j, i))
      temperature_data(j, i) = temperature_data(j, i) + 300.05
    enddo
  enddo

  !Write a restart file.
  !Open a restart file and initialize the file object.
  call open_check(open_file(fileobj, "land_compressed.nc", "overwrite", &
                            pelist=pelist, nc_format="netcdf4", is_restart=.true.))

  !Add some global attributes to the restart file.
  call register_global_attribute(fileobj, "globalatt1", real(7., kind=real64))
  call register_global_attribute(fileobj, "globalatt2", real(4., kind=real32))
  call register_global_attribute(fileobj, "globalatt3", int(3, kind=int32))
  call register_global_attribute(fileobj, "globalatt4", int(2, kind=int64))
  call register_global_attribute(fileobj, "globalatt5", "some text", str_len=len_trim("some text"))

  !Add dimensions and corresponding variables to the file.
  !cdim (compressed dimension).
  call register_axis(fileobj, "cdim", npes_corner=npes_corner, npes_nelems=npes_nelems)
  call register_field(fileobj, "cdim", "double", (/"cdim"/))
  allocate(cdim_data(nx))
  do i = 1, nx
    cdim_data(i) = 4.178 + 3.9567*i
  enddo
  call write_data(fileobj, "cdim", cdim_data)
  call register_variable_attribute(fileobj, "cdim", "units", "compressed value", str_len=len_trim("compressed value"))

  !Height.
  call register_axis(fileobj, "lev", nz)
  call register_field(fileobj, "lev", "double", (/"lev"/))
  allocate(lev_data(nz))
  do i = 1, nz
    lev_data(i) = 0.1111 + 20.236*i
  enddo
  call write_data(fileobj, "lev", lev_data)
  call register_variable_attribute(fileobj, "lev", "units", "mb", str_len=len_trim("mb"))

  !Add a restart variable.
  call register_restart_field(fileobj, "temperature", temperature_data, (/"cdim","lev "/))
  call register_variable_attribute(fileobj, "temperature", "units", "K", str_len=len_trim("K"))

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
  call write_restart(fileobj)

  !Close the file and destroy the file object.
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Prepare to read the data back in.
  deallocate(cdim_data)
  deallocate(lev_data)
  deallocate(temperature_data)

  !Read restart with an explicit axis register call.
  !Re-open the restart file and re-initialize the file object.
  call open_check(open_file(fileobj, "land_compressed.nc", "read", pelist=pelist, &
                            nc_format="netcdf4", is_restart=.true.))

  !Read the cdim variable back in and make sure it is the
  !correct size.
  ndims = get_variable_num_dimensions(fileobj, "cdim", broadcast=.true.)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, "cdim", dim_sizes, .true.)
  if (dim_sizes(1) .ne. total_nx) then
    call mpp_error(fatal, "compressed axis is not the correct size.")
  endif
  deallocate(dim_sizes)

  !Allocate a buffer and associate it with the
  !temperature variable.
  ndims = get_variable_num_dimensions(fileobj, "temperature", broadcast=.true.)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, "temperature", dim_sizes, broadcast=.true.)
  allocate(temperature_data(dim_sizes(1), dim_sizes(2)))
  deallocate(dim_sizes)
  temperature_data(:, :) = 0.

  !Read in the restart variable and close the file.
  call read_data(fileobj, "temperature", temperature_data)
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Check to see if the read-in data matches the data that was
  !written out.
  s = npes_corner(my_rank+1)
  in_chksum = mpp_chksum(temperature_data(s:s+nx-1,:), pelist=(/mpp_pe()/))
  call chksum_match(out_chksum, in_chksum, "temperature", debug)
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
  if (mpp_pe() .eq. 0) then
    write(output_unit,'(a/)') trim(green)//"passed."//trim(color_end)
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)
end subroutine land_compressed_restart_file
