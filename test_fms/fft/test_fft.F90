program test_fft
use fft_mod
integer, parameter :: lot = 2
real   , allocatable :: ain(:,:), aout(:,:)
complex, allocatable :: four(:,:)
integer :: i, j, m, n
integer :: ntrans(2) = (/ 60, 90 /)

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
    enddo
    enddo

    call fft_end
    deallocate (ain,aout,four)
  enddo

end program test_fft
