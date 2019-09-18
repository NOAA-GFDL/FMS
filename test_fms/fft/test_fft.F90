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

program test_fft
use fft_mod
use mpp_mod,       only : mpp_init, mpp_exit, mpp_error, FATAL, stdout

integer, parameter :: lot = 2
real   , allocatable :: ain(:,:), aout(:,:)
complex, allocatable :: four(:,:)
integer :: i, j, m, n
integer :: ntrans(2) = (/ 60, 90 /)

call mpp_init

! test multiple transform lengths
  do m = 1,2

  ! set up input data
    n = ntrans(m)
    allocate (ain(n+1,lot),aout(n+1,lot),four(n/2+1,lot))
    call random_number (ain(1:n,:))
    aout(1:n,:) = ain(1:n,:)

    call fft_init (n)

  ! transform grid to fourier and back
    four = fft_grid_to_fourier (aout)
    aout = fft_fourier_to_grid (four)

  ! print original and transformed
    do j=1,lot
    do i=1,n
      write (*,'(2i4,3(2x,f15.9))') j, i, ain(i,j), aout(i,j), aout(i,j)-ain(i,j)
      if (abs(aout(i,j)-ain(i,j)) > 1.0e-8) call mpp_error(FATAL, 'test_fft: the difference between the transformed and original grid is too large')

    enddo
    enddo

    call fft_end
    deallocate (ain,aout,four)
  enddo
  call mpp_exit
end program test_fft
