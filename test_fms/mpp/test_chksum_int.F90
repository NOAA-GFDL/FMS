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
!> @author Ryan Mulhall
!> @email gfdl.climate.model.info@noaa.gov
!> @brief Test mpp_chksum with mixed precision integers
!> @description Tests mpp_chksum with 8 and 4 byte integer arrays with
!> normal and distributed checksums
program test_chksum_int

  use platform_mod
  use mpp_mod, only : mpp_init, mpp_exit, mpp_pe, mpp_npes, mpp_root_pe, stdout
  use mpp_mod, only : mpp_set_stack_size, mpp_sync, mpp_init_test_init_true_only
  use mpp_mod, only : mpp_transmit, mpp_chksum, ALL_PES
  use mpp_mod, only : mpp_error, FATAL, mpp_sync_self, NOTE
  use mpp_io_mod, only: mpp_io_init

  implicit none

  integer                        :: pe, npes, root, out_unit, ierr
  integer(i8_kind), allocatable  :: data8(:), distData(:),temp(:)
  integer(i8_kind)               :: res4, res8, resDist
  integer(i4_kind), allocatable  :: data4(:)
  real, allocatable              :: rands(:)
  integer                        :: i, length

  call mpp_init(mpp_init_test_init_true_only)
  call mpp_io_init()
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
  out_unit = stdout()

  !> generate random arrays
  length = 1024
  allocate(rands(length), data8(length), data4(length), distData(length))
  call random_number(rands)
  do i = 1, length
    data8(i) = rands(i) * huge(data4(1))
    data4(i) = rands(i) * huge(data4(1))
    distData(i) = rands(i) * huge(distData(1))
  end do
  !>test mixed precision int checksums
  res4 = mpp_chksum(data4)
  res8 = mpp_chksum(data8)
  if(res4.NE.res8) then
    call mpp_error(FATAL, 'Test mpp_chksum_int: mixed precision checksums do not match')
  else
    call mpp_error(NOTE, 'Test mpp_chksum_int: mixed precision checksums match')
  endif
  !>test distributed int checksums
  call mpp_sync()
  call mpp_transmit( put_data=distData(1), plen=length, to_pe=ALL_PES, &
                     get_data=distData(1),glen=length, from_pe=root)
  call mpp_sync_self()
  allocate(temp(length/npes))
  temp = distData( pe*(length/npes)+1 : (pe+1)*(length/npes))!> distribute data for pelist
  resDist = mpp_chksum(distData(1:length), (/pe/))
  if(resDist.NE.mpp_chksum(temp)) then
    call mpp_error(FATAL, 'Test mpp_chksum_int: distributed checksums do not match')
  else
    call mpp_error(NOTE, 'Test mpp_chksum_int: distributed checksums match')
  endif
  deallocate(rands, data8, data4, distData, temp)

  call MPI_FINALIZE(ierr)

end program test_chksum_int
