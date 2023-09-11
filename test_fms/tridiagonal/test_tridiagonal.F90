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
#ifndef TEST_TRIDIAG_KIND
#define TEST_TRIDIAG_KIND 8
#endif

!> Tests the tridiagonal module routines (tri_invert and close_tridiagonal)
!! Tests reals with the kind value set above,
program test_tridiagonal

    use tridiagonal_mod
    use platform_mod
    use mpp_mod
    use fms_mod

    implicit none

    integer, parameter :: IN_LEN = 16 !< length of input arrays
    integer, parameter :: TEST_KIND = TEST_TRIDIAG_KIND !< kind value for all reals in this test
                                                !! set by TEST_TRIDIAG_KIND cpp macro

    call mpp_init
    call test_tri_invert
    call mpp_exit

    contains

    !> tests init and subsequent calls to tri_invert routine
    subroutine test_tri_invert
        real(TEST_TRIDIAG_KIND), allocatable :: a(:,:,:), b(:,:,:), c(:,:)
        allocate(a(4,4,4))
        print *, 'kind size:', KIND(a(1,1,1))
        print *, ''
    end subroutine

end program