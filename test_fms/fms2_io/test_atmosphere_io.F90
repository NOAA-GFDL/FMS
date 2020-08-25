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

program test_atmosphere_io
#ifndef use_mpp_io
use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64, error_unit, output_unit
use mpi
use fms2_io_mod
use mpp_domains_mod
use mpp_mod
use setup


type(Params) :: test_params
type(domain2d) :: domain
type(domain2d) :: domain2
integer :: err
type(domain2d), pointer :: io_domain
integer :: isc
integer :: iec
integer :: nx
integer :: jsc
integer :: jec
integer :: ny
integer :: isd
integer :: nxd
integer :: jsd
integer :: nyd
integer :: nx_east
integer :: isc_east
integer :: iec_east
integer :: ny_north
integer :: jsc_north
integer :: jec_north
type(FmsNetcdfDomainFile_t) :: fileobj
type(FmsNetcdfDomainFile_t) :: fileobjv
integer, parameter :: ntiles = 6
integer, parameter :: nt = 11
real(kind=real64), dimension(:), allocatable :: double_buffer
real(kind=real64), dimension(:,:), allocatable :: double_buffer2d
real(kind=real64), dimension(:,:,:), allocatable :: var5
real(kind=real64), dimension(:,:,:), allocatable :: var6
real(kind=real64), dimension(:,:,:), allocatable :: var7
real(kind=real64), dimension(:,:,:), allocatable :: var8
real(kind=real64), dimension(:,:,:), allocatable :: var5p
real(kind=real64), dimension(:,:,:), allocatable :: var6p
real(kind=real64), dimension(:,:,:), allocatable :: var7p
real(kind=real64), dimension(:,:,:), allocatable :: var8p
real(kind=real64), dimension(:), allocatable :: var9
integer(kind=int32), dimension(:,:), allocatable :: var10
real(kind=real64), dimension(:,:,:), allocatable :: var11
integer(kind=int64) :: var5_chksum
integer(kind=int64) :: var6_chksum
integer(kind=int64) :: var7_chksum
integer(kind=int64) :: var8_chksum
integer(kind=int64) :: var9_chksum
integer(kind=int64) :: var10_chksum
integer(kind=int64) :: var11_chksum
integer(kind=int64) :: chksum
character(len=256), dimension(:), allocatable :: string_buffer
character(len=256), dimension(:), allocatable :: string_buffer2
integer :: i
integer :: ndims
integer, dimension(:), allocatable :: dim_sizes
character(len=256) :: att
character(len=6), dimension(4) :: names
character(len=8) :: timestamp

!Initialize.
call init(test_params, ntiles)
call create_cubed_sphere_domain(test_params, domain, (/1, 1/))
call create_cubed_sphere_domain(test_params, domain2, (/1, mpp_npes()/ntiles/))
call mpi_barrier(mpi_comm_world, err)
call mpi_check(err)
call random_seed()
call fms2_io_init()

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

!Open a restart file and initialize the file object.
call open_check(open_file(fileobj, "atmosphere_io.nc", "overwrite", &
                          domain, nc_format="64bit", is_restart=.true.))
call open_check(open_virtual_file(fileobjv, domain, "atm.nc"))

!Add some global attributes to the restart file.
call register_global_attribute(fileobj, "globalatt1", real(7., kind=real64))
call register_global_attribute(fileobj, "globalatt2", real(4., kind=real32))
call register_global_attribute(fileobj, "globalatt3", int(3, kind=int32))
call register_global_attribute(fileobj, "globalatt4", int(2, kind=int64))
call register_global_attribute(fileobj, "globalatt5", "some text", str_len=9)
call register_global_attribute(fileobjv, "globalatt1", real(7., kind=real64))
call register_global_attribute(fileobjv, "globalatt2", real(4., kind=real32))
call register_global_attribute(fileobjv, "globalatt3", int(3, kind=int32))
call register_global_attribute(fileobjv, "globalatt4", int(2, kind=int64))
call register_global_attribute(fileobjv, "globalatt5", "some text", str_len=9)

!Add dimensions and corresponding variables to the file.
!Longitude (domain "x" dimension with center position).
names(1) = "lon"
call register_axis(fileobj, "lon", "x")
call register_field(fileobj, "lon", "double", names(1:1))
call create_data(double_buffer, nx)
call write_data(fileobj, "lon", double_buffer)
call register_variable_attribute(fileobj, "lon", "units", "degrees", str_len=7)
call register_axis(fileobjv, "lon", "x")
call register_field(fileobjv, "lon", "double", names(1:1))
call write_data(fileobjv, "lon", double_buffer)
call register_variable_attribute(fileobjv, "lon", "units", "degrees", str_len=7)

!Longitude2 (domain "x" dimension with east position).
names(1) = "lon2"
call register_axis(fileobj, "lon2", "x", domain_position=east)
call register_field(fileobj, "lon2", "double", names(1:1))
call create_data(double_buffer, nx+1)
call write_data(fileobj, "lon2", double_buffer)
call register_variable_attribute(fileobj, "lon2", "units", "radians", str_len=7)
call register_axis(fileobjv, "lon2", "x", domain_position=east)
call register_field(fileobjv, "lon2", "double", names(1:1))
call write_data(fileobjv, "lon2", double_buffer)
call register_variable_attribute(fileobjv, "lon2", "units", "radians", str_len=7)

!Latitude (domain "y" dimension with center position).
names(1) = "lat"
call register_axis(fileobj, "lat", "y", domain_position=center)
call register_field(fileobj, "lat", "double", names(1:1))
call create_data(double_buffer, ny)
call write_data(fileobj, "lat", double_buffer)
call register_variable_attribute(fileobj, "lat", "units", "degrees", str_len=7)
call register_axis(fileobjv, "lat", "y", domain_position=center)
call register_field(fileobjv, "lat", "double", names(1:1))
call write_data(fileobjv, "lat", double_buffer)
call register_variable_attribute(fileobjv, "lat", "units", "degrees", str_len=7)

!Latitude2 (domain "y" dimension wiht north position).
names(1) = "lat2"
call register_axis(fileobj, "lat2", "y", domain_position=north)
call register_field(fileobj, "lat2", "double", names(1:1))
call create_data(double_buffer, ny+1)
call write_data(fileobj, "lat2", double_buffer)
call register_variable_attribute(fileobj, "lat", "units2", "radians", str_len=7)
call register_axis(fileobjv, "lat2", "y", domain_position=north)
call register_field(fileobjv, "lat2", "double", names(1:1))
call write_data(fileobjv, "lat2", double_buffer)
call register_variable_attribute(fileobjv, "lat", "units2", "radians", str_len=7)

!Height.
names(1) = "lev"
call register_axis(fileobj, "lev", test_params%nz)
call register_field(fileobj, "lev", "double", names(1:1))
call create_data(double_buffer, test_params%nz)
call write_data(fileobj, "lev", double_buffer)
call register_variable_attribute(fileobj, "lev", "units", "mb", str_len=2)
call register_axis(fileobjv, "lev", test_params%nz)
call register_field(fileobjv, "lev", "double", names(1:1))
call write_data(fileobjv, "lev", double_buffer)
call register_variable_attribute(fileobjv, "lev", "units", "mb", str_len=2)


!Height2.
names(1) = "lay"
call register_axis(fileobj, "lay", test_params%nz-1)
call register_field(fileobj, "lay", "double", names(1:1))
call create_data(double_buffer, test_params%nz-1)
call write_data(fileobj, "lay", double_buffer)
call register_variable_attribute(fileobj, "lay", "units", "mb", str_len=2)
call register_axis(fileobjv, "lay", test_params%nz-1)
call register_field(fileobjv, "lay", "double", names(1:1))
call write_data(fileobjv, "lay", double_buffer)
call register_variable_attribute(fileobjv, "lay", "units", "mb", str_len=2)

!Time.
names(1) = "time"
call register_axis(fileobj, "time", unlimited)
call register_field(fileobj, "time", "float", names(1:1))
call register_variable_attribute(fileobj, "time", "units", "years", str_len=5)
call register_axis(fileobjv, "time", unlimited)
call register_field(fileobjv, "time", "float", names(1:1))
call register_variable_attribute(fileobjv, "time", "units", "years", str_len=5)

!String length.
call register_axis(fileobj, "strlen", 256)
!call register_axis(fileobjv, "strlen", 256)

!Add a non-domain-decomposed variable whose type does not match the type of the user
!defined buffer.
names(1) = "lev"
call create_data(double_buffer, test_params%nz)
call register_field(fileobj, "var1", "float", names(1:1))
call write_data(fileobj, "var1", double_buffer)
!call register_field(fileobjv, "var1", "float", names(1:1))
!call write_data(fileobjv, "var1", double_buffer)

!Add scalar and 1D string variables.
names(1) = "strlen"
names(2) = "lev"
call register_field(fileobj, "var2", "char", names(1:1))
call write_data(fileobj, "var2", "file1.nc")
call register_field(fileobj, "var3", "char", names(1:2))
allocate(string_buffer(test_params%nz))
do i = 1, test_params%nz
  string_buffer(i) = ""
  write(string_buffer(i), *) i, "string"
  string_buffer(i) = trim(adjustl(string_buffer(i)))
enddo
call write_data(fileobj, "var3", string_buffer)
!call register_field(fileobjv, "var2", "char", names(1:1))
!call write_data(fileobjv, "var2", "file1.nc")
!call register_field(fileobjv, "var3", "char", names(1:2))
!call write_data(fileobjv, "var3", string_buffer)

!Add a domain decomposed variable.
names(1) = "lon"
names(2) = "lat"
call create_data(double_buffer2d, (/nx, ny/))
call register_field(fileobj, "var4", "double", names(1:2))
call register_variable_attribute(fileobj, "var4", "units", "K", str_len=1)
call write_data(fileobj, "var4", double_buffer2d)
!call register_field(fileobjv, "var4", "double", names(1:2))
!call register_variable_attribute(fileobjv, "var4", "units", "K", str_len=1)
!call write_data(fileobjv, "var4", double_buffer2d)

!Add a domain-decomposed restart variable with center position.
names(1) = "lon"
names(2) = "lat"
names(3) = "lev"
names(4) = "time"
call create_data(var5, (/nx, ny, test_params%nz/))
call register_restart_field(fileobj, "var5", var5, names(1:4))
call register_variable_attribute(fileobj, "var5", "units", "K", str_len=1)
call create_data(var5p, (/nx, ny, test_params%nz/))
call register_restart_field(fileobjv, "var5", var5p, names(1:4))
call register_variable_attribute(fileobjv, "var5", "units", "K", str_len=1)

!Add a domain-decomposed restart variable with east position.
names(1) = "lon2"
names(2) = "lat"
names(3) = "lev"
names(4) = "time"
call create_data(var6, (/nx+1, ny, test_params%nz/))
call register_restart_field(fileobj, "var6", var6, names(1:4))
call register_variable_attribute(fileobj, "var6", "units", "K", str_len=1)
call create_data(var6p, (/nx+1, ny, test_params%nz/))
call register_restart_field(fileobjv, "var6", var6p, names(1:4))
call register_variable_attribute(fileobjv, "var6", "units", "K", str_len=1)

!Add a domain-decomposed restart variable with north position.
names(1) = "lon"
names(2) = "lat2"
names(3) = "lev"
names(4) = "time"
call create_data(var7, (/nx, ny+1, test_params%nz/))
call register_restart_field(fileobj, "var7", var7, names(1:4))
call register_variable_attribute(fileobj, "var7", "units", "K", str_len=1)
call create_data(var7p, (/nx, ny+1, test_params%nz/))
call register_restart_field(fileobjv, "var7", var7p, names(1:4))
call register_variable_attribute(fileobjv, "var7", "units", "K", str_len=1)

!Add a domain-decomposed restart variable with corner (east+north) position.
names(1) = "lon2"
names(2) = "lat2"
names(3) = "lev"
names(4) = "time"
call create_data(var8, (/nx+1, ny+1, test_params%nz/))
call register_restart_field(fileobj, "var8", var8, names(1:4))
call register_variable_attribute(fileobj, "var8", "units", "K", str_len=1)
call create_data(var8p, (/nx+1, ny+1, test_params%nz/))
call register_restart_field(fileobjv, "var8", var8p, names(1:4))
call register_variable_attribute(fileobjv, "var8", "units", "K", str_len=1)

!Add a non-domain-decomposed variable.
names(1) = "lay"
names(2) = "time"
call create_data(var9, test_params%nz-1)
call register_field(fileobj, "var9", "double", names(1:2))
call register_variable_attribute(fileobj, "var9", "units", "K", str_len=1)

!Add a domain decomposed variable.
names(1) = "lon"
names(2) = "lat"
names(3) = "time"
call create_data(var10, (/nx, ny/))
call register_field(fileobj, "var10", "int", names(1:3))
call register_variable_attribute(fileobj, "var10", "units", "K", str_len=1)

!Add a variable whose user-level buffers include space for halos.
names(1) = "lon"
names(2) = "lat"
names(3) = "lev"
names(4) = "time"
call create_data(var11, (/nxd, nyd, test_params%nz/))
call register_restart_field(fileobj, "var11", var11, names(1:4))

!Perform a "simulation" and write restart data.
do i = 1, nt
  call write_data(fileobj, "time", real(i, kind=real32), unlim_dim_level=i)
  call create_data(var5)
  call create_data(var6)
  call create_data(var7)
  call create_data(var8)
  call create_data(var9)
  call write_data(fileobj, "var9", var9, unlim_dim_level=i)
  call create_data(var10)
  call write_data(fileobj, "var10", var10, unlim_dim_level=i)
  call create_data(var11)
  call write_restart(fileobj, unlim_dim_level=i)
  call create_data(var5p)
  call create_data(var6p)
  call create_data(var7p)
  call create_data(var8p)
enddo
timestamp = "00002"
call write_new_restart(fileobjv, timestamp=timestamp)

!Store checksums for non-domain-decomposed/non-restart variables since they are not
!currently written to the output file.
!if (fileobj%is_root) then
  var9_chksum = mpp_chksum(var9, pelist=(/mpp_pe()/))
!endif
var5_chksum = mpp_chksum(var5, pelist=(/mpp_pe()/))
var6_chksum = mpp_chksum(var6, pelist=(/mpp_pe()/))
var7_chksum = mpp_chksum(var7, pelist=(/mpp_pe()/))
var8_chksum = mpp_chksum(var8, pelist=(/mpp_pe()/))
var10_chksum = mpp_chksum(var10, pelist=(/mpp_pe()/))
!var11_chksum = mpp_chksum(var11(isc-isd+1:isc-isd+1+nx, jsc-jsd+1:jsc-jsd+1+ny, :), pelist=(/mpp_pe()/))

!Close the file.
call close_file(fileobj)
call mpi_barrier(mpi_comm_world, err)
call mpi_check(err)

deallocate(double_buffer)
deallocate(double_buffer2d)
deallocate(var5)
deallocate(var6)
deallocate(var7)
deallocate(var8)
deallocate(var9)
deallocate(var10)
deallocate(var11)

!Check if a non-existent file exists (just to see if this feature
!works.
if (open_file(fileobj, "atmosphere.foobar.nc", "read", domain, &
              nc_format="64bit", is_restart=.true.)) then
  call mpp_error(FATAL, "Found non-existent file.")
endif

!Re-open the restart file and re-initialize the file object.
call open_check(open_file(fileobj, "atmosphere_io.nc", "read", domain2, &
                          nc_format="64bit", is_restart=.true.))

!Get the sizes of the I/O compute and data domains.
io_domain => mpp_get_io_domain(domain2)
if (.not. associated(io_domain)) then
  call mpp_error(fatal, "I/O domain is not associated.")
endif
call mpp_get_compute_domain(io_domain, xbegin=isc, xend=iec, xsize=nx, ybegin=jsc, yend=jec, ysize=ny)
call mpp_get_data_domain(io_domain, xbegin=isd, xsize=nxd, ybegin=jsd, ysize=nyd)
call mpp_get_compute_domain(io_domain, xbegin=isc_east, xend=iec_east, xsize=nx_east, position=east)
call mpp_get_compute_domain(io_domain, ybegin=jsc_north, yend=jec_north, ysize=ny_north, position=north)

!Associate the lon and lat dimensions with the
!"x" and "y" dimensions of the domain.
call register_axis(fileobj, "lon", "x", domain_position=center)
call register_axis(fileobj, "lat", "y")
call register_axis(fileobj, "lon2", "x", domain_position=east)
call register_axis(fileobj, "lat2", "y", domain_position=north)

!Read in and check the string variables.
call get_dimension_size(fileobj, "lev", ndims)
allocate(string_buffer2(ndims))
call read_data(fileobj, "var3", string_buffer2)
do i = 1, ndims
  if (trim(string_buffer2(i)) .ne. trim(string_buffer(i))) then
    call mpp_error(fatal, "Did not read in strings correctly.")
  endif
enddo
call read_data(fileobj, "var3", string_buffer2(1), corner=ndims)
if (trim(string_buffer2(1)) .ne. trim(string_buffer(ndims))) then
  call mpp_error(fatal, "Scalar from array of strings read failed.")
endif
call read_data(fileobj, "var2", string_buffer2(ndims))
if (trim(string_buffer2(ndims)) .ne. "file1.nc") then
  call mpp_error(fatal, "Did not read in filename variable correctly.")
endif
deallocate(string_buffer)
deallocate(string_buffer2)

!Read in and compare chksums of non-restart variables.
ndims = get_variable_num_dimensions(fileobj, "var9", broadcast=.true.)
allocate(dim_sizes(ndims))
call get_variable_size(fileobj, "var9", dim_sizes, broadcast=.true.)
allocate(var9(dim_sizes(1)))
call read_data(fileobj, "var9", var9, unlim_dim_level=nt)
!if (fileobj%is_root) then
!  chksum = mpp_chksum(var9, pelist=(/mpp_pe()/))
!  if (chksum .ne. var9_chksum) then
!    call mpp_error(fatal, "checksum for var9 does not match.")
!  endif
!endif
deallocate(var9)
deallocate(dim_sizes)

!Add a domain decomposed variable.
ndims = get_variable_num_dimensions(fileobj, "var10", broadcast=.true.)
allocate(dim_sizes(ndims))
call get_variable_size(fileobj, "var10", dim_sizes, broadcast=.true.)
allocate(var10(nx, ny))
call read_data(fileobj, "var10", var10, unlim_dim_level=nt)
chksum = mpp_chksum(var10, pelist=(/mpp_pe()/))
if (chksum .ne. var10_chksum) then
  call mpp_error(fatal, "checksum for var10 does not match.")
endif
deallocate(var10)
deallocate(dim_sizes)

!Allocate buffers and register restart variables.
ndims = get_variable_num_dimensions(fileobj, "var5", broadcast=.true.)
allocate(dim_sizes(ndims))
call get_variable_size(fileobj, "var5", dim_sizes, broadcast=.true.)
allocate(var5(nx, ny, dim_sizes(3)))
call register_restart_field(fileobj, "var5", var5)
deallocate(dim_sizes)

ndims = get_variable_num_dimensions(fileobj, "var6", broadcast=.true.)
allocate(dim_sizes(ndims))
call get_variable_size(fileobj, "var6", dim_sizes, broadcast=.true.)
allocate(var6(nx+1, ny, dim_sizes(3)))
call register_restart_field(fileobj, "var6", var6)
deallocate(dim_sizes)

ndims = get_variable_num_dimensions(fileobj, "var7", broadcast=.true.)
allocate(dim_sizes(ndims))
call get_variable_size(fileobj, "var7", dim_sizes, broadcast=.true.)
allocate(var7(nx, ny+1, dim_sizes(3)))
call register_restart_field(fileobj, "var7", var7)
deallocate(dim_sizes)

ndims = get_variable_num_dimensions(fileobj, "var8", broadcast=.true.)
allocate(dim_sizes(ndims))
call get_variable_size(fileobj, "var8", dim_sizes, broadcast=.true.)
allocate(var8(nx+1, ny+1, dim_sizes(3)))
call register_restart_field(fileobj, "var8", var8)
deallocate(dim_sizes)

ndims = get_variable_num_dimensions(fileobj, "var11", broadcast=.true.)
allocate(dim_sizes(ndims))
call get_variable_size(fileobj, "var11", dim_sizes, broadcast=.true.)
allocate(var11(nxd, nyd, dim_sizes(3)))
call register_restart_field(fileobj, "var11", var11)
deallocate(dim_sizes)

!Read in the restart data.
call read_restart(fileobj, unlim_dim_level=nt)


chksum = mpp_chksum(var5, pelist=(/mpp_pe()/))
if (chksum .ne. var5_chksum) then
  call mpp_error(fatal, "checksum for var 5 does not match.")
else
  call mpp_error(warning, "checksum for var 5 does match.")
endif

chksum = mpp_chksum(var6, pelist=(/mpp_pe()/))
if (chksum .ne. var6_chksum) then
  call mpp_error(fatal, "checksum for var 6 does not match.")
else
  call mpp_error(warning, "checksum for var 6 does match.")
endif

chksum = mpp_chksum(var7, pelist=(/mpp_pe()/))
if (chksum .ne. var7_chksum) then
  call mpp_error(fatal, "checksum for var 7 does not match.")
else
  call mpp_error(warning, "checksum for var 7 does match.")
endif

chksum = mpp_chksum(var8, pelist=(/mpp_pe()/))
if (chksum .ne. var8_chksum) then
  call mpp_error(fatal, "checksum for var 8 does not match.")
else
  call mpp_error(warning, "checksum for var 8 does match.")
endif





var5p = 0.
var6p = 0.
var7p = 0.
var8p = 0.
call read_new_restart(fileobjv, timestamp=timestamp)

!Close the file.
call close_file(fileobj)
call close_file(fileobjv)
call mpi_barrier(mpi_comm_world, err)
call mpi_check(err)

!Clean up.
if (mpp_pe() .eq. 0) then
  write(output_unit,'(a/)') &
    "Atmosphere (6-tile domain decomposed) file test ... " &
    //trim(green)//"passed."//trim(color_end)
endif
call mpi_barrier(mpi_comm_world, err)
call mpi_check(err)
deallocate(var5)
deallocate(var6)
deallocate(var7)
deallocate(var8)
deallocate(var11)
deallocate(var5p)
deallocate(var6p)
deallocate(var7p)
deallocate(var8p)
call cleanup(test_params)




































#ifdef FOOBAR

!Read in variable attribute.
call get_variable_attribute(fileobj, "temperature", "units", att)
if (trim(att) .ne. "K") then
  call mpp_error(fatal, "failed to read in string variable attribute.")
endif

#endif






#endif
end program test_atmosphere_io
