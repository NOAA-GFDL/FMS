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
!! Arguments A, B, and C are optional, and are saved as module variables
!! if one recalls tri_invert without changing (A,B,C)

!> @file
!> @brief File for @ref tridiagonal_mod

!> @addtogroup tridiagonal_mod
!> @{
module tridiagonal_mod

!--------------------------------------------------------------------------
real,    private, allocatable, dimension(:,:,:) :: e,g,cc
real,    private, allocatable, dimension(:,:)   :: bb
logical, private :: init_tridiagonal = .false.
!--------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------

!> @brief Sets up and solves the tridiagonal system of equations
!!
!> For simplicity, A and C are assumed to be dimensioned the same size
!! as B, D, and X, although any input values for A(N) and C(1) are ignored.
!! There are no checks to make sure the sizes agree.
!!
!! The value of A(N) is modified on output, and B and C are unchanged.
subroutine tri_invert(x,d,a,b,c)

implicit none

real, intent(out), dimension(:,:,:) :: x !< Solution to the tridiagonal system of equations
real, intent(in),  dimension(:,:,:) :: d !< The right-hand side term, see the schematic above.
real, optional,    dimension(:,:,:) :: a,b,c !< Left hand side terms(see schematic above).
                                             !! If not provided, values from last call are used

real, dimension(size(x,1),size(x,2),size(x,3)) :: f
integer :: k

if(present(a)) then
  init_tridiagonal = .true.

  if(allocated(e))     deallocate(e)
  if(allocated(g))     deallocate(g)
  if(allocated(bb))    deallocate(bb)
  if(allocated(cc))    deallocate(cc)
  allocate(e (size(x,1),size(x,2),size(x,3)))
  allocate(g (size(x,1),size(x,2),size(x,3)))
  allocate(bb(size(x,1),size(x,2)))
  allocate(cc(size(x,1),size(x,2),size(x,3)))

  e(:,:,1) = - a(:,:,1)/b(:,:,1)
  a(:,:,size(x,3)) = 0.0

  do  k= 2,size(x,3)
    g(:,:,k) = 1.0/(b(:,:,k)+c(:,:,k)*e(:,:,k-1))
    e(:,:,k) = - a(:,:,k)*g(:,:,k)
  end do
  cc = c
  bb = 1.0/b(:,:,1)

end if

! if(.not.init_tridiagonal) error

f(:,:,1) =  d(:,:,1)*bb
do k= 2, size(x,3)
  f(:,:,k) = (d(:,:,k) - cc(:,:,k)*f(:,:,k-1))*g(:,:,k)
end do

x(:,:,size(x,3)) = f(:,:,size(x,3))
do k = size(x,3)-1,1,-1
  x(:,:,k) = e(:,:,k)*x(:,:,k+1)+f(:,:,k)
end do

return
end subroutine tri_invert

!-----------------------------------------------------------------

!> @brief Releases memory used by the solver
subroutine close_tridiagonal

implicit none

deallocate(e)
deallocate(g)
deallocate(bb)
deallocate(cc)

return
end subroutine close_tridiagonal

!----------------------------------------------------------------

end module tridiagonal_mod

! <INFO>

!   <BUG>
!     Optional arguments A,B,C have no intent declaration,
!     so the default intent is inout. The value of A(N) is modified
!     on output, and B and C are unchanged.
!   </BUG>
!   <NOTE>
!       The following private allocatable arrays save the relevant information
!  if one recalls tri_invert without changing (A,B,C):
!  <PRE>
!        allocate ( e  (size(x,1), size(x,2), size(x,3)) )
!        allocate ( g  (size(x,1), size(x,2), size(x,3)) )
!        allocate ( cc (size(x,1), size(x,2), size(x,3)) )
!        allocate ( bb (size(x,1), size(x,2)) )
! </PRE>
!  This storage is deallocated when close_tridiagonal is called.
!   </NOTE>
!   <FUTURE>
!     Maybe a cleaner version?
!   </FUTURE>

! </INFO>
!> @}
! close documentation grouping
