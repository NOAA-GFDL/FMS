
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

program test_drifters_input
  use drifters_input_mod
  use fms_mod, only : fms_init, fms_end
  use mpp_mod, only : mpp_error, FATAL, stdout

  implicit none
  character(len=128) :: ermesg
  integer :: i

  type(drifters_input_type) :: obj

  call fms_init()

  call drifters_input_new(obj, 'test_drifters_input.nc', ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)

  print *,'field_names:'
  do i = 1, size(obj%field_names)
     print *,trim(obj%field_names(i))
  enddo

  print *,'velocity_names:'
  do i = 1, size(obj%velocity_names)
     print *,trim(obj%velocity_names(i))
  enddo

  print *,'ids = ', obj%ids

  print *,'positions: '
  do i = 1, size(obj%positions, 2)
     print *,obj%positions(:,i)
  enddo

  call drifters_input_del(obj, ermesg)

  call fms_end()
end program test_drifters_input
