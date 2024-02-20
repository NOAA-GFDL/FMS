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

!> @brief  Checks the output for when running with a field that has a mask that changes
!! over time
program check_var_masks
  use fms_mod,           only: fms_init, fms_end
  use mpp_mod
  use fms2_io_mod

  implicit none

  type(FmsNetcdfFile_t)              :: fileobj
  integer                            :: ntimes
  integer                            :: nx
  integer                            :: ny
  real,   allocatable                :: vardata(:,:)
  real                               :: ans_var_mask
  real                               :: ans_var
  integer                            :: i, j

  call fms_init()

  if (.not. open_file(fileobj, "test_var_masks.nc", "read")) &
    call mpp_error(FATAL, "unable to open test_var_masks.nc for reading")

  call get_dimension_size(fileobj, "time", ntimes)
  if (ntimes .ne. 1) call mpp_error(FATAL, "time is not the correct size!")

  call get_dimension_size(fileobj, "x", nx)
  if (nx .ne. 360) call mpp_error(FATAL, "x is not the correct size!")

  call get_dimension_size(fileobj, "y", ny)
  if (ny .ne. 180) call mpp_error(FATAL, "y is not the correct size!")

  allocate(vardata(nx,ny))

  ans_var_mask = 0.
  ans_var = 0.
  call read_data(fileobj, "ua", vardata)
  do i = 1, 24
    ans_var = ans_var + real(i)
    if (mod(i,2) .ne. 0) ans_var_mask = ans_var_mask + real(i)
  enddo
  ans_var = ans_var / 24
  ans_var_mask = ans_var_mask / 12

  do i = 1, nx
    do j = 1, ny
      if (i .eq. 1 .and. j .eq. 1) then
        if (vardata(i,j) .ne. ans_var_mask) &
          call mpp_error(FATAL, "ua is not the expected result for the masked point")
      else
        if (vardata(i,j) .ne. ans_var) &
          call mpp_error(FATAL, "ua is not the expected result")
      endif
    enddo
  enddo

  call close_file(fileobj)
  call fms_end()
end program check_var_masks
