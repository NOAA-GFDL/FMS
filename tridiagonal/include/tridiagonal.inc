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

!> @addtogroup tridiagonal_mod
!> @{

!> @brief Sets up and solves the tridiagonal system of equations
!!
!> For simplicity, A and C are assumed to be dimensioned the same size
!! as B, D, and X, although any input values for A(N) and C(1) are ignored.
!! There are no checks to make sure the sizes agree.
!!
!! The value of A(N) is modified on output, and B and C are unchanged.
!! 
!! For mixed precision, this routine uses the kind size macro(FMS_TRID_KIND_) to determine
!! which module variables are used/stored. This means a,b, and c values will only be stored for calls
!! of the same real kind value. 
subroutine TRI_INVERT_(x,d,a,b,c)

    real(FMS_TRID_KIND_), intent(out), dimension(:,:,:) :: x !< Solution to the tridiagonal system of equations
    real(FMS_TRID_KIND_), intent(in),  dimension(:,:,:) :: d !< The right-hand side term, see the schematic above.
    real(FMS_TRID_KIND_), optional,    dimension(:,:,:) :: a,b,c !< Left hand side terms(see schematic above).
                                                !! If not provided, values from last call are used

    real(FMS_TRID_KIND_), dimension(size(x,1),size(x,2),size(x,3)) :: f

    integer :: k

    if(present(a)) then
        !$OMP SINGLE
        init_tridiagonal = .true.
        if(allocated(TRID_REAL_TYPE%e))     deallocate(TRID_REAL_TYPE%e)
        if(allocated(TRID_REAL_TYPE%g))     deallocate(TRID_REAL_TYPE%g)
        if(allocated(TRID_REAL_TYPE%bb))    deallocate(TRID_REAL_TYPE%bb)
        if(allocated(TRID_REAL_TYPE%cc))    deallocate(TRID_REAL_TYPE%cc)
        if(allocated(TRID_REAL_TYPE))       deallocate(TRID_REAL_TYPE)
        allocate(TRID_REAL_TYPE)
        allocate(TRID_REAL_TYPE%e (size(x,1),size(x,2),size(x,3)))
        allocate(TRID_REAL_TYPE%g (size(x,1),size(x,2),size(x,3)))
        allocate(TRID_REAL_TYPE%bb(size(x,1),size(x,2)))
        allocate(TRID_REAL_TYPE%cc(size(x,1),size(x,2),size(x,3)))
        !$OMP END SINGLE

        TRID_REAL_TYPE%e(:,:,1) = - a(:,:,1) / b(:,:,1)
        a(:,:,size(x,3)) = real(0.0, FMS_TRID_KIND_)

        do  k= 2,size(x,3)
            TRID_REAL_TYPE%g(:,:,k) = real(1.0, FMS_TRID_KIND_)/(b(:,:,k)+c(:,:,k)*TRID_REAL_TYPE%e(:,:,k-1))
            TRID_REAL_TYPE%e(:,:,k) = - a(:,:,k)* TRID_REAL_TYPE%g(:,:,k)
        end do
        TRID_REAL_TYPE%cc = c
        TRID_REAL_TYPE%bb = 1.0/b(:,:,1)

    end if

    if(.not.init_tridiagonal) call mpp_error(FATAL, 'tri_invert: a,b,and c args not provided or previously calculated.')

    f(:,:,1) =  d(:,:,1)*TRID_REAL_TYPE%bb
    do k= 2, size(x,3)
        f(:,:,k) = (d(:,:,k) - TRID_REAL_TYPE%cc(:,:,k)*f(:,:,k-1))*TRID_REAL_TYPE%g(:,:,k)
    end do

    x(:,:,size(x,3)) = f(:,:,size(x,3))
    do k = size(x,3)-1,1,-1
        x(:,:,k) = TRID_REAL_TYPE%e(:,:,k)*x(:,:,k+1)+f(:,:,k)
    end do

    return
end subroutine TRI_INVERT_