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
!> @file
!! @brief unit test for mpp_broadcast
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests the subroutines
!! mpp_broadcast_real4_2d to mpp_broadcast_real4_5d,
!! mpp_broadcast_real8_2d to mpp_broadcast_real8_5d,
!! mpp_broadcast_i4_2d to mpp_broadcast_i4_5d,
!! mpp_broadcast_i8_2d to mpp_broadcast_i8_5d, and mpp_broadcast_char

program test_mpp_broadcast

  use platform_mod
  use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod, only : mpp_error, mpp_broadcast, FATAL

  implicit none

  integer :: ierr !< Used by MPI_FINALIZE

  call mpp_init(test_level=mpp_init_test_peset_allocated)

  !> tests mpp_broadcast_*D_I4
  call test_broadcast_2D_I4()
  call test_broadcast_3D_I4()
  call test_broadcast_4D_I4()
  call test_broadcast_5D_I4()
  !> tests mpp_broadcast_*D_I8
  call test_broadcast_2D_I8()
  call test_broadcast_3D_I8()
  call test_broadcast_4D_I8()
  call test_broadcast_5D_I8()
  !> tests mpp_broadcast_*D_R4
  call test_broadcast_2D_R4()
  call test_broadcast_3D_R4()
  call test_broadcast_4D_R4()
  call test_broadcast_5D_R4()
  !> tests mpp_broadcast_*D_R8
  call test_broadcast_2D_R8()
  call test_broadcast_3D_R8()
  call test_broadcast_4D_R8()
  call test_broadcast_5D_R8()
  !> tests mpp_broadcast_char
  call test_broadcast_char()

  call MPI_FINALIZE(ierr)

contains
!>
!> test mpp_broadcast_2d_i4
!>
subroutine test_broadcast_2D_I4()

  implicit none

  integer, parameter :: NN = 3
  integer(i4_kind), parameter  :: zero = 0, one=1

  integer :: n, m
  integer(i4_kind) :: p, r(NN,NN), k(NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        p = p + one
        k(m,n) = p
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           if(r(m,n) == k(m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
        enddo
     enddo
 endif

 call mpp_broadcast(r, NN*NN, mpp_root_pe())

 !--- after broadcast, r and k should be the same
 do n=1, NN
    do m=1, NN
       if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
    enddo
 enddo

end subroutine test_broadcast_2D_I4
!>
!> test mpp_broadcast_3d_i4
!>
subroutine test_broadcast_3D_i4()

  implicit none

  integer, parameter :: NN = 3
  integer(i4_kind), parameter  :: zero = 0, one=1

  integer :: i, n, m
  integer(i4_kind) :: p, r(NN,NN,NN), k(NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           p = p + one
           k(i,m,n) = p
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) == k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
        enddo
     enddo
  enddo

end subroutine test_broadcast_3D_i4
!>
!> test mpp_broadcoast_4D_i4
!>
subroutine test_broadcast_4D_i4()

  implicit none

  integer, parameter :: NN = 3
  integer(i4_kind), parameter  :: zero = 0, one=1

  integer :: i, j, n, m
  integer(i4_kind) :: p, r(NN,NN,NN,NN), k(NN,NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              p = p + one
              k(j,i,m,n) = p
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) == k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                   call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_4D_i4
!>
!> test mpp_broadcast_5d_i4
!>
subroutine test_broadcast_5D_I4()

  implicit none

  integer, parameter :: NN = 3
  integer(i4_kind), parameter  :: zero = 0, one=1

  integer :: i, j, l, n, m
  integer(i4_kind) :: p, r(NN,NN,NN,NN,NN), k(NN,NN,NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 p = p + one
                 k(l,j,i,m,n) = p
              enddo
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) == k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_5D_I4
!>
!> test_broadcast_2d_i8
!>
subroutine test_broadcast_2D_I8()

  implicit none

  integer, parameter :: NN = 3
  integer(i8_kind), parameter  :: zero = 0, one=1

  integer :: n, m
  integer(i8_kind) :: p, r(NN,NN), k(NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        p = p + one
        k(m,n) = p
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           if(r(m,n) == k(m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
        enddo
     enddo
 endif

 call mpp_broadcast(r, NN*NN, mpp_root_pe())

 !--- after broadcast, r and k should be the same
 do n=1, NN
    do m=1, NN
       if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
    enddo
 enddo

end subroutine test_broadcast_2D_I8
!>
!> test_broadcast_3D_i8
!>
subroutine test_broadcast_3D_i8()

  implicit none

  integer, parameter :: NN = 3
  integer(i8_kind), parameter  :: zero = 0, one=1

  integer :: i, n, m
  integer(i8_kind) :: p, r(NN,NN,NN), k(NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           p = p + one
           k(i,m,n) = p
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) == k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
        enddo
     enddo
  enddo

end subroutine test_broadcast_3D_i8
!>
!> test mpp_broadcast_4d_i8
!>
subroutine test_broadcast_4D_i8()

  implicit none

  integer, parameter :: NN = 3
  integer(i8_kind), parameter  :: zero = 0, one=1

  integer :: i, j, n, m
  integer(i8_kind) :: p, r(NN,NN,NN,NN), k(NN,NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              p = p + one
              k(j,i,m,n) = p
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) == k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                   call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_4D_i8
!>
!> test mpp_broadcast_5d_i8
!>
subroutine test_broadcast_5D_I8()

  implicit none

  integer, parameter :: NN = 3
  integer(i8_kind), parameter  :: zero = 0, one=1

  integer :: i, j, l, n, m
  integer(i8_kind) :: p, r(NN,NN,NN,NN,NN), k(NN,NN,NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 p = p + one
                 k(l,j,i,m,n) = p
              enddo
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) == k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_5D_I8
!>
!> test mpp_broadcast_2d_r4
!>
subroutine test_broadcast_2D_R4()

  implicit none

  integer, parameter :: NN = 3
  real(r4_kind), parameter  :: zero = 0., one=1.

  integer :: n, m
  real(r4_kind) :: p, r(NN,NN), k(NN,NN)

  p=zero
  do n = 1, NN
     do m = 1, NN
        p = p + one
        k(m,n) = p
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n = 1, NN
        do m = 1, NN
           if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
        enddo
     enddo
  else
     do n = 1, NN
        do m = 1, NN
           if(r(m,n) == k(m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
        enddo
     enddo
 endif

 call mpp_broadcast(r, NN*NN, mpp_root_pe())

 !--- after broadcast, r and k should be the same
 do n = 1, NN
    do m = 1, NN
       if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
    enddo
 enddo

end subroutine test_broadcast_2D_R4
!>
!> test mpp_broadcast_3d_r4
!>
subroutine test_broadcast_3D_R4()

  implicit none

  integer, parameter :: NN = 3
  real(r4_kind), parameter  :: zero = 0., one=1.

  integer :: i, n, m
  real(r4_kind) :: p, r(NN,NN,NN), k(NN,NN,NN)

  p=zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           p = p + one
           k(i,m,n) = p
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) == k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
        enddo
     enddo
  enddo

end subroutine test_broadcast_3D_R4


!> test mpp_broadcoast_4D_R4
subroutine test_broadcast_4D_R4()

  implicit none

  integer, parameter :: NN = 3
  real(r4_kind), parameter  :: zero = 0., one=1.

  integer :: i, j, n, m
  real(r4_kind) :: p, r(NN,NN,NN,NN), k(NN,NN,NN,NN)

  p=zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              p = p + one
              k(j,i,m,n) = p
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) == k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                   call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_4D_R4

!> test mpp_broadcast_5d_r4
subroutine test_broadcast_5D_R4()

  implicit none

  integer, parameter :: NN = 3
  real(r4_kind), parameter  :: zero = 0., one=1.

  integer :: i, j, l, n, m
  real(r4_kind) :: p, r(NN,NN,NN,NN,NN), k(NN,NN,NN,NN,NN)

  p=zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 p = p + one
                 k(l,j,i,m,n) = p
              enddo
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) == k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_5D_R4


!> test mpp_broadcast_2d_r8
subroutine test_broadcast_2D_R8()

  implicit none

  integer, parameter :: NN = 3
  real(r8_kind), parameter  :: zero = 0., one=1.

  integer :: n, m
  real(r8_kind) :: p, r(NN,NN), k(NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        p = p + one
        k(m,n) = p
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
        enddo
     enddo
  else
     do n = 1, NN
        do m = 1, NN
           if(r(m,n) == k(m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        if(r(m,n) .NE. k(m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
     enddo
  enddo

end subroutine test_broadcast_2D_R8
!>
!> test mpp_broadcast_3d_r8
!>
subroutine test_broadcast_3D_R8()

  implicit none

  integer, parameter :: NN = 3
  real(r8_kind), parameter  :: zero = 0., one=1.

  integer :: i, n, m
  real(r8_kind) :: p, r(NN,NN,NN), k(NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           p = p + one
           k(i,m,n) = p
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              if(r(i,m,n) == k(i,m,n)) call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           if(r(i,m,n) .NE. k(i,m,n)) call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
        enddo
     enddo
  enddo

end subroutine test_broadcast_3D_R8
!>
!> test mpp_broadcast_4d_R8
!>
subroutine test_broadcast_4D_R8()
  implicit none

  integer, parameter :: NN = 3
  real(r8_kind), parameter  :: zero = 0., one=1.

  integer :: i, j, n, m
  real(r8_kind) :: p, r(NN,NN,NN,NN), k(NN,NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              p = p + one
              k(j,i,m,n) = p
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 if(r(j,i,m,n) == k(j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              if(r(j,i,m,n) .NE. k(j,i,m,n)) &
                   call mpp_error(FATAL, "test_broadcast: after broadcast, r should equal k")
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_4D_R8
!>
!> test mpp_broadcast_5d_r8
!>
subroutine test_broadcast_5D_R8()

  implicit none

  integer, parameter :: NN = 3
  real(r8_kind), parameter  :: zero = 0., one=1.

  integer :: i, j, l, n, m
  real(r8_kind) :: p, r(NN,NN,NN,NN,NN), k(NN,NN,NN,NN,NN)

  p = zero
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 p = p + one
                 k(l,j,i,m,n) = p
              enddo
           enddo
        enddo
     enddo
  enddo

  r = k
  if(mpp_pe() .NE. mpp_root_pe()) r = zero

  !--- comparing array r and k. r and k are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on root_pe, r should equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  else
     do n=1, NN
        do m=1, NN
           do i=1, NN
              do j=1, NN
                 do l=1, NN
                    if(r(l,j,i,m,n) == k(l,j,i,m,n)) &
                         call mpp_error(FATAL, "test_broadcast: on non root_pes, r should not equal k")
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  call mpp_broadcast(r, NN*NN*NN*NN*NN, mpp_root_pe())

  !--- after broadcast, r and k should be the same
  do n=1, NN
     do m=1, NN
        do i=1, NN
           do j=1, NN
              do l=1, NN
                 if(r(l,j,i,m,n) .NE. k(l,j,i,m,n)) &
                      call mpp_error(FATAL, "test_broadcast: after broadcast, r should equalk")
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine test_broadcast_5D_R8
!>
!> test mpp_broadcast_char
!>
subroutine test_broadcast_char()

  implicit none

  integer, parameter :: NN = 3
  integer, parameter :: STRINGSIZE = 256
  character(len=STRINGSIZE), dimension(NN) :: textA, textB
  integer :: n

  textA(1) = "This is line 1 "
  textA(2) = "Here comes the line 2 "
  textA(3) = "Finally is line 3 "
  do n = 1, NN
     textB(n) = TextA(n)
  enddo

  if(mpp_pe() .NE. mpp_root_pe()) then
     do n =1, NN
        textA(n) = ""
     enddo
  endif

  !--- comparing textA and textB. textA and textB are supposed to be
  !different on pe other than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n = 1, NN
        if(textA(n) .NE. textB(n)) call mpp_error(FATAL, "test_broadcast: on root_pe, textA should equal textB")
     enddo
  else
     do n = 1, NN
        if(textA(n) == textB(n)) call mpp_error(FATAL, "test_broadcast: on root_pe, textA should not equal textB")
     enddo
  endif
  call mpp_broadcast(textA, STRINGSIZE, mpp_root_pe())
  !--- after broadcast, textA and textB should be the same
  do n = 1, NN
     if(textA(n) .NE. textB(n)) call mpp_error(FATAL, "test_broadcast: after broadcast, textA should equal textB")
  enddo

end subroutine test_broadcast_char

end program test_mpp_broadcast
