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
!> @defgroup tridiagonal_mod tridiagonal_mod
!> @ingroup tridiagonal
!> @brief Solves a tridiagonal system of equations.
!!
!> The following schematic represents the system of equations solved,
!! where X is the solution.
!! <PRE>
!!     | B(1)  A(1)   0     0                .......            0    |  |X(1)|   |D(1)|
!!     | C(2)  B(2)  A(2)   0                .......            0    |  |X(2)|   |D(2)|
!!     |  0    C(3)  B(3)  A(3)  0           .......            0    |  | .. |   | .. |
!!     |  ..........................................                 |  | .. | = | .. |
!!     |  ..........................................                 |  | .. |   | .. |
!!     |                                  C(N-2) B(N-2) A(N-2)  0    |  | .. |   | .. |
!!     |                                    0    C(N-1) B(N-1) A(N-1)|  | .. |   | .. |
!!     |                                    0      0    C(N)   B(N)  |  |X(N)|   |D(N)|
!!
!! </PRE>
!!  To solve this system
!! <PRE>
!!   call tri_invert(X,D,A,B,C)
!!
!!       real, intent(out), dimension(:,:,:) :: X
!!       real, intent(in),  dimension(:,:,:) :: D
!!       real, optional,    dimension(:,:,:) :: A,B,C
!! </PRE>
!! For simplicity (?), A and C are assumed to be dimensioned the same size
!! as B, D, and X, although any input values for A(N) and C(1) are ignored.
!! (some checks are needed here)
!!
!! If A is not present, it is assumed that the matrix (A,B.C) has not been changed
!! since the last call to tri_invert.
!!
!! To release memory,
!! <PRE>
!!    call close_tridiagonal
!! </PRE>
!!
!!
!! Arguments A, B, and C are optional, and are saved as module variables
!! if one recalls tri_invert without changing (A,B,C)
!!
!! @note
!!     Optional arguments A,B,C have no intent declaration,
!!     so the default intent is inout. The value of A(N) is modified
!!     on output, and B and C are unchanged.
!!
!!  The following private allocatable arrays save the relevant information
!!  if one recalls tri_invert without changing (A,B,C):
!!  <PRE>
!!        allocate ( e  (size(x,1), size(x,2), size(x,3)) )
!!        allocate ( g  (size(x,1), size(x,2), size(x,3)) )
!!        allocate ( cc (size(x,1), size(x,2), size(x,3)) )
!!        allocate ( bb (size(x,1), size(x,2)) )
!! </PRE>
!!  This storage is deallocated when close_tridiagonal is called.

!> @addtogroup tridiagonal_mod
!> @{
module tridiagonal_mod

    use platform_mod, only: r4_kind, r8_kind
    use mpp_mod,      only: mpp_error, FATAL
    implicit none

    type :: tridiag_reals_r4
        real(r4_kind), private, allocatable, dimension(:,:,:) :: e, g, cc
        real(r4_kind), private, allocatable, dimension(:,:)   :: bb
    end type

    type :: tridiag_reals_r8
        real(r8_kind), private, allocatable, dimension(:,:,:) :: e, g, cc
        real(r8_kind), private, allocatable, dimension(:,:)   :: bb
    end type

    type(tridiag_reals_r4) :: tridiag_r4
    type(tridiag_reals_r8) :: tridiag_r8

    !! allocated when a,b,c are passed to tri_invert
    logical, private :: init_tridiagonal_r4 = .false.
    logical, private :: init_tridiagonal_r8 = .false.

    !> Interface to solve tridiagonal systems of equations for either kind value.
    !! Since this relies on the state of module variables (unless A,B,C are specified)
    !! the values stored are distinct for each kind call unless the added optional argument store_both_kinds
    !! is true
    interface tri_invert
        module procedure tri_invert_r4
        module procedure tri_invert_r8
    end interface

    public :: tri_invert

    contains

    !> @brief Releases memory used by the solver
    subroutine close_tridiagonal
        !$OMP SINGLE
        if(.not. init_tridiagonal_r4 .and. .not. init_tridiagonal_r8) return
        if(allocated(tridiag_r4%e)) deallocate(tridiag_r4%e)
        if(allocated(tridiag_r4%g)) deallocate(tridiag_r4%g)
        if(allocated(tridiag_r4%cc)) deallocate(tridiag_r4%cc)
        if(allocated(tridiag_r4%bb)) deallocate(tridiag_r4%bb)
        if(allocated(tridiag_r8%e)) deallocate(tridiag_r8%e)
        if(allocated(tridiag_r8%g)) deallocate(tridiag_r8%g)
        if(allocated(tridiag_r8%cc)) deallocate(tridiag_r8%cc)
        if(allocated(tridiag_r8%bb)) deallocate(tridiag_r8%bb)
        init_tridiagonal_r4 = .false.; init_tridiagonal_r8 = .false.
        !$OMP END SINGLE
        return
    end subroutine close_tridiagonal

#include "tridiagonal_r4.fh"
#include "tridiagonal_r8.fh"

end module tridiagonal_mod

!> @}
! close documentation grouping