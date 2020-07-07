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

program test_mpp_broadcast

  use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod, only : mpp_error, mpp_broadcast, FATAL

  integer :: ierr !< Used by MPI_FINALIZE

  call mpp_init(test_level=mpp_init_test_peset_allocated)

  call test_broadcast_2D()
  call test_broadcast_char()

  call MPI_FINALIZE(ierr)
contains 

subroutine test_broadcast_2D()
  integer, parameter :: ARRAYSIZE = 3
  integer :: n, m, p
  real :: r(3,3), k(3,3)

  p=0;
  do n = 1, ARRAYSIZE
    do m = 1, ARRAYSIZE
       p = p + 1
       k(n, m) = p
       r(n, m) = k(n, m)
    enddo
  enddo

  if(mpp_pe() .NE. mpp_root_pe()) then
    do n =1, ARRAYSIZE
       r(:, n) = 0
    enddo
  endif

  !--- comparing array m and n. m and n are supposed to be different on pe other
  !than root_pe
  if(mpp_pe() == mpp_root_pe()) then
    do n = 1, ARRAYSIZE
       do m = 1, ARRAYSIZE
          if(r(n, m) .NE. k(n, m)) call mpp_error(FATAL, "test_broadcast: on root_pe, m should equal n")
       enddo
    enddo
  else
    do n = 1, ARRAYSIZE
       do m = 1, ARRAYSIZE
          if(r(n, m) == k(n, m)) call mpp_error(FATAL, "test_broadcast: on non root_pes, m should equal n")
       enddo
    enddo
  endif

  call mpp_broadcast(r, ARRAYSIZE*ARRAYSIZE, mpp_root_pe())

  !--- after broadcast, m and n should be the same
  do n = 1, ARRAYSIZE
     do m =1, ARRAYSIZE
        if(r(n, m) .NE. k(n, m)) call mpp_error(FATAL, "test_broadcast: after broadcast, m should equal n")
     enddo
  enddo

end subroutine test_broadcast_2D

subroutine test_broadcast_char()
  integer, parameter :: ARRAYSIZE = 3
  integer, parameter :: STRINGSIZE = 256
  character(len=STRINGSIZE), dimension(ARRAYSIZE) :: textA, textB
  integer :: n

  textA(1) = "This is line 1 "
  textA(2) = "Here comes the line 2 "
  textA(3) = "Finally is line 3 "
  do n = 1, ARRAYSIZE
     textB(n) = TextA(n)
  enddo

  if(mpp_pe() .NE. mpp_root_pe()) then
     do n =1, ARRAYSIZE
        textA(n) = ""
     enddo
  endif

  !--- comparing textA and textB. textA and textB are supposed to be
  !different on pe other than root_pe
  if(mpp_pe() == mpp_root_pe()) then
     do n = 1, ARRAYSIZE
        if(textA(n) .NE. textB(n)) call mpp_error(FATAL, "test_broadcast: on root_pe, textA should equal textB")
     enddo
  else
     do n = 1, ARRAYSIZE
        if(textA(n) == textB(n)) call mpp_error(FATAL, "test_broadcast: on root_pe, textA should not equal textB")
     enddo
  endif
  call mpp_broadcast(textA, STRINGSIZE, mpp_root_pe())
  !--- after broadcast, textA and textB should be the same
  do n = 1, ARRAYSIZE
     if(textA(n) .NE. textB(n)) call mpp_error(FATAL, "test_broadcast: after broadcast, textA should equal textB")
  enddo

end subroutine test_broadcast_char

end program test_mpp_broadcast

