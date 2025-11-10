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

program test_chunksizes
  use fms2_io_mod,  only: open_file, close_file, register_axis, register_restart_field, write_restart, &
                          unlimited, fmsnetcdffile_t, read_restart
  use mpp_mod,      only: mpp_error, FATAL
  use fms_mod,      only: fms_init, fms_end
  use platform_mod, only: r8_kind

  implicit none

  integer, parameter     :: dim_len = 24                           !< The dimension length
  integer                :: i, j, k                                !< For do loops
  type(fmsnetcdffile_t)  :: fileobj                                !< FMS2_io fileobj
  character (len = 120) :: my_format(3)                            !< Array of formats to try.
  character (len = 120) :: dimnames(4)                             !< Array of dimension names
  integer, dimension(4)  :: chunksizes                             !< The chunksizes to use
  real(kind=r8_kind)     :: vardata_in(dim_len, dim_len, dim_len)  !< The data to write
  real(kind=r8_kind)     :: vardata_out(dim_len, dim_len, dim_len) !< The data read in

  call fms_init()

  my_format(1) = '64bit'
  my_format(2) = 'classic'
  my_format(3) = 'netcdf4'

  dimnames = (/"dim1", "dim2", "dim3", "dim4"/)

  chunksizes = (/dim_len, dim_len, dim_len, 1/)

  do i = 1, dim_len
    do j = 1, dim_len
      do k = 1, dim_len
        vardata_in(i, j, k) = real(k, kind=r8_kind)/13.
      enddo
    enddo
  enddo

  !< Loop through each of the file formats and write out a file
  do i = 1, 4
    if (i .ne. 4) then
      if (.not. open_file(fileobj, "test_chunksizes_"//trim(my_format(i))//".nc", "overwrite", &
        is_restart=.true., nc_format=my_format(i))) &
      call mpp_error(FATAL, "Error opening the file to write:test_chunksizes_"//trim(my_format(i))//".nc")
    else
      if (.not. open_file(fileobj, "test_chunksizes.nc", "overwrite", is_restart=.true.)) &
      call mpp_error(FATAL, "Error opening the file to write:test_chunksizes.nc")
    endif

    call register_axis(fileobj, "dim1", dim_len)
    call register_axis(fileobj, "dim2", dim_len)
    call register_axis(fileobj, "dim3", dim_len)
    call register_axis(fileobj, "dim4", unlimited)
    !< The chunksizes will ignored for non netcdf4 file formats
    call register_restart_field(fileobj, "var1", vardata_in, dimnames, chunksizes=chunksizes)
    call write_restart(fileobj)
    call close_file(fileobj)
  enddo

  !< Loop through each of the file formats and read out a file
  do i = 1, 4
    if (i .ne. 4) then
      if (.not. open_file(fileobj, "test_chunksizes_"//trim(my_format(i))//".nc", "read", &
        is_restart=.true., nc_format=my_format(i))) &
      call mpp_error(FATAL, "Error opening the file to read:test_chunksizes_"//trim(my_format(i))//".nc")
    else
      if (.not. open_file(fileobj, "test_chunksizes.nc", "read", is_restart=.true.)) &
      call mpp_error(FATAL, "Error opening the file to read:test_chunksizes.nc")
    endif

    vardata_out = -999.
    call register_restart_field(fileobj, "var1", vardata_out, dimnames, chunksizes=chunksizes)
    call read_restart(fileobj)
    call close_file(fileobj)

    if (sum(vardata_out) .ne. sum(vardata_in)) &
      call mpp_error(FATAL, "Error reading the data back from file")
  enddo

  call fms_end()

end program
