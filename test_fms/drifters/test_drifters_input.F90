
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

program test_drifters_input
#ifdef use_drifters
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
#endif
end program test_drifters_input
