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
!> @file
!! @brief unit test for mpp_alltoall
!! @author MiKyung Lee
!! @description
!! test sending/receiving 1 element so that, for example for npes=4,
!! process0: [ 0, 1, 2, 3]  --alltoall--> [0,10,20,30]
!! process1: [10,11,12,13]  --alltoall--> [1,11,21,31]
!! process2: [20,21,22,23]  --alltoall--> [2,12,22,32]
!! process3: [30,31,32,33]  --alltoall--> [3,13,23,33]
!! and test sending/receiving more than 1 element so that, for example, for npes=4 and Nsend=nrecv=2
!! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoall--> [0,1,10,11,20,21,30,31]
!! process1: [10,11,12,13,14,15,16,17]  --alltoall--> [2,3,12,13,22,23,32,33]
!! process2: [20,21,22,23,24,25,26,27]  --alltoall--> [4,5,14,15,24,25,34,35]
!! process3: [30,31,32,33,34,35,36,37]  --alltoall--> [6,7,16,17,26,27,36,37]
!! https://www.olcf.ornl.gov/wp-content/uploads/2018/06/intro_to_HPC_intro_to_mpi.pdf

program test_mpp_alltoall

  use platform_mod
  use mpp_mod, only : mpp_init, mpp_init_test_requests_allocated, mpp_init_test_peset_allocated, mpp_error, FATAL
  use mpp_mod, only : mpp_pe, mpp_npes, mpp_alltoall
  use mpp_mod, only : mpp_type_create, mpp_type, mpp_byte

  implicit none

  integer :: npes, ierr

    !> initialize MPI
    call mpp_init( test_level=mpp_init_test_requests_allocated )

    !> get total number of pe's
    npes = mpp_npes()

    !> call tests
    call test_mpp_alltoall_real4(npes)
    call test_mpp_alltoall_real8(npes)
    call test_mpp_alltoall_int4(npes)
    call test_mpp_alltoall_int8(npes)

    call test_mpp_alltoallv_real4(npes)
    call test_mpp_alltoallv_real8(npes)
    call test_mpp_alltoallv_int4(npes)
    call test_mpp_alltoallv_int8(npes)

    call test_mpp_alltoallw_real4(npes)
    call test_mpp_alltoallw_real8(npes)
    call test_mpp_alltoallw_int4(npes)
    call test_mpp_alltoallw_int8(npes)

    call MPI_FINALIZE(ierr)


  contains

  !>
  !> test mpp_alltoall for real4
  !>

  subroutine test_mpp_alltoall_real4(npes)

    implicit none

    integer, intent(in) :: npes

    real(r4_kind), parameter :: zero = 0., one = 1.

    integer :: pe, ierr, i, ii, N, isend, jsend, irecv, nsend, nrecv
    real(r4_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving up to npes elements.  can set up to 9 elements.
    nsend = npes  ;  nrecv = nsend
    if ( npes > 9 ) then
       nsend = 9  ;  nrecv = 9
    end if

    do isend=1, nsend

       !> allocate sbuf (senddata), rbuf (receivedata)
       N = isend*npes - 1
       allocate( sbuf(0:N), rbuf(0:N) )

       !> initialize receiving array
       rbuf = -one

       !> intialize sending array
       do i=0, N
          sbuf(i) = real( 10*pe+i, kind=r4_kind )
       end do

       !> number of elements to send and receive
       irecv = isend

       !> call mpp_alltoall to send/receive one element
       call mpp_alltoall( sbuf, isend, rbuf, irecv )

       !> check
       ii = 0
       do i=0, (npes-1)
          do jsend=0, isend-1
             if( rbuf(ii) .ne. real( 10*i+isend*pe+jsend, kind=r4_kind ) ) then
                write(*,'("PE #",i3,"element",i4,"Expected",f6.0,"but received",f6.0)') pe, ii, real(10*i+nsend*pe+jsend), rbuf(ii)
                call mpp_error(FATAL, 'test_mpp_alltoall failed')
             end if
             ii = ii + 1
          end do
       end do

       deallocate( sbuf, rbuf )

    end do

  end subroutine test_mpp_alltoall_real4

  !>
  !> test mpp_alltoall for real8
  !>

  subroutine test_mpp_alltoall_real8(npes)

    implicit none

    integer, intent(in) :: npes

    real(r8_kind), parameter :: zero = 0., one=1.

    integer :: pe, ierr, i, ii, N, isend, jsend, irecv, nsend, nrecv
    real(r8_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving up to npes elements.  can set up to 9 elements.
    nsend = npes  ;  nrecv = nsend
    if ( npes > 9 ) then
       nsend = 9  ;  nrecv = 9
    end if

    do isend=1, nsend

       !> allocate sbuf (senddata), rbuf (receivedata)
       N = isend*npes - 1
       allocate( sbuf(0:N), rbuf(0:N) )

       !> initialize receiving array
       rbuf = - one

       !> intialize sending array
       do i=0, N
          sbuf(i) = real( 10*pe+i, kind=r8_kind )
       end do

       !> number of elements to send and receive
       irecv = isend

       !> call mpp_alltoall to send/receive one element
       call mpp_alltoall( sbuf, isend, rbuf, irecv )

       !> check
       ii = 0
       do i=0, (npes-1)
          do jsend=0, isend-1
             if( rbuf(ii) .ne. real( 10*i+isend*pe+jsend, kind=r8_kind ) ) then
                write(*,'("PE #",i3,"element",i4,"Expected",f6.0,"but received",f6.0)') pe, ii, real(10*i+nsend*pe+jsend), rbuf(ii)
                call mpp_error(FATAL, 'test_mpp_alltoall failed')
             end if
             ii = ii + 1
          end do
       end do

       deallocate( sbuf, rbuf )

    end do

  end subroutine test_mpp_alltoall_real8

  !>
  !> test mpp_alltoall for int4
  !>

  subroutine test_mpp_alltoall_int4(npes)

    implicit none

    integer, intent(in) :: npes

    integer(i4_kind), parameter :: zero = 0, one=1

    integer :: pe, ierr, i, ii, N, isend, jsend, irecv, nsend, nrecv
    integer(i4_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving up to npes elements.  can set up to 9 elements.
    nsend = npes  ;  nrecv = nsend
    if ( npes > 9 ) then
       nsend = 9  ;  nrecv = 9
    end if

    do isend=1, nsend

       !> allocate sbuf (senddata), rbuf (receivedata)
       N = isend*npes - 1
       allocate( sbuf(0:N), rbuf(0:N) )

       !> initialize receiving array
       rbuf = - one

       !> intialize sending array
       do i=0, N
          sbuf(i) = real( 10*pe+i, kind=i4_kind )
       end do

       !> number of elements to send and receive
       irecv = isend

       !> call mpp_alltoall to send/receive one element
       call mpp_alltoall( sbuf, isend, rbuf, irecv )

       !> check
       ii = 0
       do i=0, (npes-1)
          do jsend=0, isend-1
             if( rbuf(ii) .ne. real( 10*i+isend*pe+jsend, kind=i4_kind ) ) then
                write(*,'("PE #",i3,"element",i4,"Expected",i6,"but received",i6)') pe, ii, real(10*i+nsend*pe+jsend), rbuf(ii)
                call mpp_error(FATAL, 'test_mpp_alltoall failed')
             end if
             ii = ii + 1
          end do
       end do

       deallocate( sbuf, rbuf )

    end do

  end subroutine test_mpp_alltoall_int4

  !>
  !> test mpp_alltoall for int8
  !>

  subroutine test_mpp_alltoall_int8(npes)

    implicit none

    integer, intent(in) :: npes

    integer(i8_kind), parameter :: zero = 0, one=1

    integer :: pe, ierr, i, ii, N, isend, jsend, irecv, nsend, nrecv
    integer(i8_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving up to npes elements.  can set up to 9 elements.
    nsend = npes  ;  nrecv = nsend
    if ( npes > 9 ) then
       nsend = 9  ;  nrecv = 9
    end if

    do isend=1, nsend

       !> allocate sbuf (senddata), rbuf (receivedata)
       N = isend*npes - 1
       allocate( sbuf(0:N), rbuf(0:N) )

       !> initialize receiving array
       rbuf = - one

       !> intialize sending array
       do i=0, N
          sbuf(i) = real( 10*pe+i, kind=i8_kind )
       end do

       !> number of elements to send and receive
       irecv = isend

       !> call mpp_alltoall to send/receive one element
       call mpp_alltoall( sbuf, isend, rbuf, irecv )

       !> check
       ii = 0
       do i=0, (npes-1)
          do jsend=0, isend-1
             if( rbuf(ii) .ne. real( 10*i+isend*pe+jsend, kind=i8_kind ) ) then
                write(*,'("PE #",i3,"element",i4,"Expected",i6,"but received",i6)') pe, ii, real(10*i+nsend*pe+jsend), rbuf(ii)
                call mpp_error(FATAL, 'test_mpp_alltoall failed')
             end if
             ii = ii + 1
          end do
       end do

       deallocate( sbuf, rbuf )

    end do

  end subroutine test_mpp_alltoall_int8

  !>
  !> test mpp_alltoallv for real4
  !>

  subroutine test_mpp_alltoallv_real4(npes)

    implicit none

    integer, intent(in) :: npes

    real(r4_kind) :: zero = 0., one = 1.

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    real(r4_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N = npes - 1

    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    ssize = 1  ;  rsize = 1

    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    do i=0, N
       sbuf(i) = real( 10*pe+i, kind=r4_kind )
    end do
    rbuf = -one

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl )

    !> check
    do i=0, N
       if ( rbuf(i).ne.real(10*i+pe, kind=r4_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

    !>send one element, receive one
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    ssize = 1  ;  rsize = 1

    deallocate( sbuf, rbuf )
    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    do i=0, N
       sbuf(2*i)   = real( 10*pe+2*i, kind=r4_kind )
       sbuf(2*i+1) = real( 10*pe+2*i+1, kind=r4_kind )
    end do

    rbuf = real(-1.0, kind=r4_kind )

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.real(10*i+2*pe, kind=r4_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
       if ( rbuf(2*i+1).ne.-one ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

  end subroutine test_mpp_alltoallv_real4

  !>
  !> test mpp_alltoallv for real8
  !>

  subroutine test_mpp_alltoallv_real8(npes)

    implicit none

    integer, intent(in) :: npes

    real(r8_kind) :: zero = 0., one = 1.

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    real(r8_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N = npes - 1

    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    ssize = 1  ;  rsize = 1

    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    do i=0, N
       sbuf(i) = real( 10*pe+i, kind=r8_kind )
    end do
    rbuf = -one

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl )

    !> check
    do i=0, N
       if ( rbuf(i).ne.real(10*i+pe, kind=r8_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

    !>send one element, receive one
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    ssize = 1  ;  rsize = 1

    deallocate( sbuf, rbuf )
    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    do i=0, N
       sbuf(2*i)   = real( 10*pe+2*i, kind=r8_kind )
       sbuf(2*i+1) = real( 10*pe+2*i+1, kind=r8_kind )
    end do

    rbuf = real(-1.0, kind=r8_kind )

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.real(10*i+2*pe, kind=r8_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
       if ( rbuf(2*i+1).ne.-one ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

  end subroutine test_mpp_alltoallv_real8

  !>
  !> test mpp_alltoallv for int4
  !>

  subroutine test_mpp_alltoallv_int4(npes)

    implicit none

    integer, intent(in) :: npes

    integer(i4_kind) :: zero = 0, one = 1

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    real(i4_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N = npes - 1

    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    ssize = 1  ;  rsize = 1

    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    do i=0, N
       sbuf(i) = real( 10*pe+i, kind=i4_kind )
    end do
    rbuf = -one

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl )

    !> check
    do i=0, N
       if ( rbuf(i).ne.real(10*i+pe, kind=i4_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

    !>send one element, receive one
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    ssize = 1  ;  rsize = 1

    deallocate( sbuf, rbuf )
    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    do i=0, N
       sbuf(2*i)   = real( 10*pe+2*i, kind=i4_kind )
       sbuf(2*i+1) = real( 10*pe+2*i+1, kind=i4_kind )
    end do

    rbuf = -one

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.real(10*i+2*pe, kind=i4_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
       if ( rbuf(2*i+1).ne.-one ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

  end subroutine test_mpp_alltoallv_int4

  !>
  !> test mpp_alltoallv for int4
  !>

  subroutine test_mpp_alltoallv_int8(npes)

    implicit none

    integer, intent(in) :: npes

    integer(i8_kind) :: zero = 0, one = 1

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    real(i8_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N = npes - 1

    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    ssize = 1  ;  rsize = 1

    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    do i=0, N
       sbuf(i) = real( 10*pe+i, kind=i8_kind )
    end do
    rbuf = -one

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl )

    !> check
    do i=0, N
       if ( rbuf(i).ne.real(10*i+pe, kind=i8_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

    !>send one element, receive one
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    ssize = 1  ;  rsize = 1

    deallocate( sbuf, rbuf )
    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    do i=0, N
       sbuf(2*i)   = real( 10*pe+2*i, kind=i8_kind )
       sbuf(2*i+1) = real( 10*pe+2*i+1, kind=i8_kind )
    end do

    rbuf = -one

    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.real(10*i+2*pe, kind=i8_kind) ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
       if ( rbuf(2*i+1).ne.-one ) call mpp_error( FATAL, 'test mpp_alltoallv fail' )
    end do

  end subroutine test_mpp_alltoallv_int8

  !>
  !> test mpp_alltoallw_real4
  !>

  subroutine test_mpp_alltoallw_real4(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n=9
    integer, parameter :: byte4 = 4
    real(r4_kind), parameter :: zero = 0. , one = 1.

    integer :: pe, i, j, jj, jjj, k, kk
    real(r4_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    real(r4_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    real(r4_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)


    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf data arrays
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = real( pe*1000 + i*100 + j*10 + k, kind=r4_kind )
          end do
       end do
    end do
    rbuf = - one

    !>
    !> test each PE sending a column of length subsize_k, and starting from sbuf(start_k,start_j,start_i)
    !>

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=real( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=r4_kind )
             if ( i.gt.subsize_i )                           answer=-one
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one
             !if( pe==1 ) write(*,*) i,j,k, rbuf(k,j,i), answer
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with columns' )
          end do
       end do
    end do

    !>
    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !>

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=real( (k-1)*1000 + i*100 + j*10 + pe+1, kind=r4_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rows' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r4_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 2 subarrays' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r4_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 3 subarrays' )
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_real4

  !>
  !> test mpp_alltoallw_real8
  !>

  subroutine test_mpp_alltoallw_real8(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n=9
    integer, parameter :: byte8 = 8
    real(r8_kind), parameter :: zero = 0. , one = 1.

    integer :: pe, i, j, jj, jjj, k, kk
    real(r8_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    real(r8_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    real(r8_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)


    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf data arrays
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = real( pe*1000 + i*100 + j*10 + k, kind=r8_kind )
          end do
       end do
    end do
    rbuf = - one

    !>
    !> test each PE sending a column of length subsize_k, and starting from sbuf(start_k,start_j,start_i)
    !>

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=real( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=r8_kind )
             if ( i.gt.subsize_i )                           answer=-one
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one
             !if( pe==1 ) write(*,*) i,j,k, rbuf(k,j,i), answer
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with columns' )
          end do
       end do
    end do

    !>
    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !>

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=real( (k-1)*1000 + i*100 + j*10 + pe+1, kind=r8_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rows' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r8_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 2 subarrays' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r8_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 3 subarrays' )
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_real8

  !>
  !> test mpp_alltoallw_int4
  !>

  subroutine test_mpp_alltoallw_int4(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n=9
    integer, parameter :: byte4 = 4
    integer(i4_kind), parameter :: zero = 0 , one = 1

    integer :: pe, i, j, jj, jjj, k, kk
    integer(i4_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    integer(i4_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    integer(i4_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)


    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf data arrays
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = int( pe*1000 + i*100 + j*10 + k, kind=i4_kind )
          end do
       end do
    end do
    rbuf = - one

    !>
    !> test each PE sending a column of length subsize_k, and starting from sbuf(start_k,start_j,start_i)
    !>

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=real( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=i4_kind )
             if ( i.gt.subsize_i )                           answer=-one
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one
             !if( pe==1 ) write(*,*) i,j,k, rbuf(k,j,i), answer
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with columns' )
          end do
       end do
    end do

    !>
    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !>

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=int( (k-1)*1000 + i*100 + j*10 + pe+1, kind=i4_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rows' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i4_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 2 subarrays' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i4_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 3 subarrays' )
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_int4

  !>
  !> test mpp_alltoallw_int8
  !>

  subroutine test_mpp_alltoallw_int8(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n=9
    integer, parameter :: byte8 = 8
    integer(i8_kind), parameter :: zero = 0 , one = 1

    integer :: pe, i, j, jj, jjj, k, kk
    integer(i8_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    integer(i8_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    integer(i8_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)


    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf data arrays
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = int( pe*1000 + i*100 + j*10 + k, kind=i8_kind )
          end do
       end do
    end do
    rbuf = - one

    !>
    !> test each PE sending a column of length subsize_k, and starting from sbuf(start_k,start_j,start_i)
    !>

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=real( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=i8_kind )
             if ( i.gt.subsize_i )                           answer=-one
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one
             !if( pe==1 ) write(*,*) i,j,k, rbuf(k,j,i), answer
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with columns' )
          end do
       end do
    end do

    !>
    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !>

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=int( (k-1)*1000 + i*100 + j*10 + pe+1, kind=i8_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rows' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i8_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 2 subarrays' )
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i8_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) call mpp_error( FATAL, 'error in MPP_ALLTOALLW with rank 3 subarrays' )
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_int8

end program test_mpp_alltoall
