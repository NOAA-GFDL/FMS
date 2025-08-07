
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

program test_drifters_core
#ifdef use_drifters

  use drifters_core_mod
  use fms_mod, only : fms_init, fms_end
  use mpp_mod, only : mpp_error, FATAL, stdout
  use platform_mod

  implicit none
  type(drifters_core_type) :: drf
  integer :: ier, nd, npdim, i, j, np
  character(128) :: ermesg
  integer :: npa
  real   , allocatable :: positions(:,:), positions_to_add(:,:)

  call fms_init()

  ! c-tor/d-tor tests
  nd    = 3
  npdim = 2
  call drifters_core_new(drf, nd, npdim, ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)
  call drifters_core_del(drf, ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)
  call drifters_core_new(drf, nd, npdim, ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)

  call drifters_core_print(drf, ermesg)

  npdim = 10
  call drifters_core_resize(drf, npdim, ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)
  call drifters_core_print(drf, ermesg)

  np = 7
  allocate(positions(nd,np))
  positions(1,:) = (/0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0/) ! x
  positions(2,:) = (/0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1/) ! y
  positions(3,:) = (/0.2, 1.2, 2.2, 3.2, 4.2, 5.2, 6.2/) ! z
  call drifters_core_set_positions(drf, positions, ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)
  call drifters_core_print(drf, ermesg)

  ! remove more particles than are added
  npa = 2
  allocate(positions_to_add(nd,npa))
  positions_to_add(1,:) = (/100.0, 200.0/)
  positions_to_add(2,:) = (/100.1, 200.1/)
  positions_to_add(3,:) = (/100.2, 200.2/)
  call drifters_core_remove_and_add(drf, (/2, 6, 1/), &
     & (/ 1001, 1002 /), &
     & positions_to_add, &
     & ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)
  call drifters_core_print(drf, ermesg)
  deallocate(positions_to_add)

  ! add more particles than are removed
  npa = 3
  allocate(positions_to_add(nd,npa))
  positions_to_add(1,:) = (/1000.0, 2000.0, 3000.0/)
  positions_to_add(2,:) = (/1000.1, 2000.1, 3000.1/)
  positions_to_add(3,:) = (/1000.2, 2000.2, 3000.2/)
  call drifters_core_remove_and_add(drf, (/3,1/), &
     & (/ 1003, 1004, 1005 /), &
     & positions_to_add,  &
     & ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)
  call drifters_core_print(drf, ermesg)
  deallocate(positions_to_add)

  ! add particles requiring resizing
  npa = 10
  allocate(positions_to_add(nd,npa))
  positions_to_add(1,:) = (/100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 10000.0/)
  positions_to_add(2,:) = (/100.1, 200.1, 300.1, 400.1, 500.1, 600.1, 700.1, 800.1, 900.1, 10000.1/)
  positions_to_add(3,:) = (/100.2, 200.2, 300.2, 400.2, 500.2, 600.2, 700.2, 800.2, 900.2, 10000.2/)
  call drifters_core_remove_and_add(drf, (/3,1,5,2/), &
     & (/ (1010+i, i=1,npa) /), &
     & positions_to_add,  &
     & ermesg)
  if(ermesg/='') call mpp_error(FATAL, ermesg)
  call drifters_core_print(drf, ermesg)
  deallocate(positions_to_add)

!!$  call test_circle(ier)
!!$  !call test_3d(ier)
!!$
!!$  if(ier/=0) then
!!$     print *,'Test unit failed ier=', ier
!!$  else
!!$     print *,'Sucessful test ier=', ier
!!$  end if
  call fms_end()
#endif
end program test_drifters_core
