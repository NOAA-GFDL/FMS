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
!> Tests the tridiagonal module routines (tri_invert and close_tridiagonal)
!! Tests reals with the kind value set above,
program test_tridiagonal

    use tridiagonal_mod
    use platform_mod
    use mpp_mod
    use fms_mod

    implicit none

    integer, parameter :: IN_LEN = 8 !< length of input arrays
    integer, parameter :: kindl = TEST_FMS_KIND_ !< kind value for all reals in this test
                                                !! set by TEST_FMS_KIND_ cpp macro
    real(TEST_FMS_KIND_), allocatable :: d(:,:,:), x(:,:,:), ref_array(:,:,:)
    real(TEST_FMS_KIND_), allocatable :: a(:,:,:), b(:,:,:), c(:,:,:)
    real(r4_kind), allocatable :: d_r4(:,:,:), x_r4(:,:,:)
    real(r4_kind), allocatable :: a_r4(:,:,:), b_r4(:,:,:), c_r4(:,:,:)
    real(r8_kind), allocatable :: d_r8(:,:,:), x_r8(:,:,:)
    real(r8_kind), allocatable :: a_r8(:,:,:), b_r8(:,:,:), c_r8(:,:,:)
    integer :: i, end, ierr, io
    real(TEST_FMS_KIND_) :: k
    ! nml
    logical :: do_error_check = .false.
    namelist / test_tridiagonal_nml/ do_error_check

    call mpp_init

    read (input_nml_file, test_tridiagonal_nml, iostat=io)
    ierr = check_nml_error (io, 'test_tridiagonal_nml')

    ! allocate input and output arrays
    allocate(d(1,1,IN_LEN))
    allocate(a(1,1,IN_LEN))
    allocate(b(1,1,IN_LEN))
    allocate(c(1,1,IN_LEN))
    allocate(x(1,1,IN_LEN))

    !! simple test with only 1 coeff
    a = 0.0_kindl
    b = 1.0_kindl
    c = 0.0_kindl
    d = 5.0_kindl
    call tri_invert(x, d, a, b, c)
    if(any(x .ne. 5.0_kindl)) call mpp_error(FATAL, "test_tridiagonal: invalid results for 1 coefficient check")
    !! check with stored data arrays
    d = -5.0_kindl
    call tri_invert(x, d)
    if(any(x .ne. -5.0_kindl)) call mpp_error(FATAL, "test_tridiagonal: invalid results for 1 coefficient check")

    ! test with a,b,c
    ! 0.5x(n-2) + x(n-1) + 0.5x(n) = 1
    !
    ! x(n) = k * [4, 1, 3, 2, 2, 3, 1, 4]
    !        k * [8 , 1, 7, 2, 6, .. ] = k *(-n/2 + ((n%2)*arr_length/2))
    a = 0.5_kindl
    b = 1.0_kindl
    c = 0.5_kindl
    d = 1.0_kindl
    call tri_invert(x, d, a, b, c)
    ! set up reference answers
    k = 1.0_kindl/(IN_LEN+1.0_kindl) * 2.0_kindl
    allocate(ref_array(1,1,IN_LEN))
    do i=1, IN_LEN/2
      end=IN_LEN-i+1
      if(mod(i, 2) .eq. 1) then
        ref_array(1,1,i) = real(-(i/2) + (mod(i,2)* IN_LEN/2), kindl)
        ref_array(1,1,end) = real(-(i/2) + (mod(i,2)* IN_LEN/2), kindl)
      else
        ref_array(1,1,i) = real(i/2, kindl)
        ref_array(1,1,end) = real(i/2, kindl)
      endif
    enddo
    ref_array = ref_array * k
    ! check
    do i=1, IN_LEN
      if(ABS(x(1,1,i) - ref_array(1,1,i)) .gt. 1.0e-6_kindl) then
        print *, i, x(1,1,i), ref_array(1,1,i)
        call mpp_error(FATAL, "test_tridiagonal: failed reference check for tri_invert")
      endif
    enddo
    !! check with stored data arrays
    d = -1.0_kindl
    ref_array = ref_array * -1.0_kindl
    call tri_invert(x, d)
    do i=1, IN_LEN
      if(ABS(x(1,1,i) - ref_array(1,1,i)) .gt. 1.0e-6_kindl) then
        print *, i, x(1,1,i), ref_array(1,1,i)
        call mpp_error(FATAL, "test_tridiagonal: failed reference check for tri_invert with saved values")
      endif
    enddo
    call close_tridiagonal()

    !! tests for module state across kinds
    !! default keeps stored values separate depending on kind
    !! store_both_kinds argument can be specified to store both r4 and r8 kinds
    if(kindl .eq. r8_kind) then
      allocate(a_r4(1,1,IN_LEN), b_r4(1,1,IN_LEN), c_r4(1,1,IN_LEN))
      allocate(d_r4(1,1,IN_LEN), x_r4(1,1,IN_LEN))
      a_r4 = 0.0_r4_kind; b_r4 = 1.0_r4_kind; c_r4 = 0.0_r4_kind
      d_r4 = 5.0_r4_kind; x_r4 = 0.0_r4_kind
      a = 0.0_kindl; b = 2.0_kindl; c = 0.0_kindl
      d = 5.0_kindl
      ! default, module variables distinct per kind
      call tri_invert(x_r4, d_r4, a_r4, b_r4, c_r4)
      ! conditionally errors here for calling with unallocated a/b/c for kind
      if( do_error_check ) call tri_invert(x, d)
      call tri_invert(x, d, a, b, c)
      ! check both values are correct from prior state
      call tri_invert(x_r4, d_r4)
      call tri_invert(x, d)
      if(any(x_r4 .ne. 5.0_r4_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r4 kind result")
      if(any(x .ne. 2.5_r8_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r8 kind result")
      call close_tridiagonal()
      ! run with storing for both kinds
      call tri_invert(x_r4, d_r4, a_r4, b_r4, c_r4, store_both_kinds=.true.)
      call tri_invert(x_r4, d_r4)
      call tri_invert(x, d)
      if(any(x_r4 .ne. 5.0_r4_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r4 kind result")
      if(any(x .ne. 5.0_r8_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r8 kind result")
    else
      allocate(a_r8(1,1,IN_LEN), b_r8(1,1,IN_LEN), c_r8(1,1,IN_LEN))
      allocate(d_r8(1,1,IN_LEN), x_r8(1,1,IN_LEN))
      a_r8 = 0.0_r8_kind; b_r8 = 1.0_r8_kind; c_r8 = 0.0_r8_kind
      d_r8 = 5.0_r8_kind; x_r8 = 0.0_r8_kind
      a = 0.0_kindl; b = 2.0_kindl; c = 0.0_kindl
      d = 5.0_kindl
      ! default, module variables distinct per kind
      call tri_invert(x_r8, d_r8, a_r8, b_r8, c_r8)
      ! conditionally errors here for calling with unallocated a/b/c for kind
      if( do_error_check ) call tri_invert(x, d)
      call tri_invert(x, d, a, b, c)
      ! check both values are correct from prior state
      call tri_invert(x_r8, d_r8)
      call tri_invert(x, d)
      if(any(x_r8 .ne. 5.0_r8_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r8 kind result")
      if(any(x .ne. 2.5_r8_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r8 kind result")
      call close_tridiagonal()
      ! run with storing for both kinds
      call tri_invert(x_r8, d_r8, a_r8, b_r8, c_r8, store_both_kinds=.true.)
      call tri_invert(x_r8, d_r8)
      call tri_invert(x, d)
      if(any(x_r8 .ne. 5.0_r8_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r8 kind result")
      if(any(x .ne. 5.0_r8_kind)) call mpp_error(FATAL, "test_tridiagonal: invalid r8 kind result")
    endif

    call close_tridiagonal()

    call mpp_exit

end program
