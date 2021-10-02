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
! nf95 -r8 -g -I ~/regression/ia64/23-Jun-2005/CM2.1U_Control-1990_E1.k32pe/include/ -D_TEST_CLOUD_INTERPOLATOR -D_F95 cloud_interpolator.F90

#define _FLATTEN(A) reshape((A), (/size((A))/) )

!> @defgroup cloud_interpolator_mod cloud_interpolator_mod
!> @ingroup drifters
!! @brief Cloud interpolation routines for use in @ref drifters_mod

!> @file
!> @brief File for @ref cloud_interpolator_mod

!> @addtogroup cloud_interpolator_mod
!> @{
MODULE cloud_interpolator_mod
  implicit none
  private

  public :: cld_ntrp_linear_cell_interp, cld_ntrp_locate_cell, cld_ntrp_get_cell_values
  public :: cld_ntrp_expand_index, cld_ntrp_contract_indices

! Include variable "version" to be written to log file.
#include<file_version.h>
real, parameter           :: tol = 10.0*epsilon(1.)

CONTAINS

!...............................................................................
!> Get expanded list of indices from contracted index
!> @param Ic contracted index
!> @param[out] ie(:) expanded list of indices
!> @param[out] ier error flag, non zero if operation unsuccessful
pure subroutine cld_ntrp_expand_index(Ic, ie, ier)
    integer, intent(in)  ::  Ic
    integer, intent(out) ::  ie(:)
    integer, intent(out) ::  ier

    integer j, nd

    ier =  0
    nd  = size(ie) ! dimension

    if(Ic >= 2**nd) then
       ie  = -1
       ier = 1 ! error
       return
    endif

    do j = 1, nd
       ie(j) = mod(Ic/2**(j-1), 2)
    end do

  end subroutine cld_ntrp_expand_index

!...............................................................................
!> Contract list of indices to an single integer
!> @param ie(:) expanded list of indices
!> @param[out] Ic contracted index
!> @param[out] ier error flag, non zero if operation unsuccessful
pure subroutine cld_ntrp_contract_indices(ie, Ic, ier)
    integer, intent(in) ::  ie(:)
    integer, intent(out)  ::  Ic
    integer, intent(out) ::  ier

    integer j, nd

    ier = 0
    nd  = size(ie) ! dimension

    Ic = ie(nd)
    do j = nd-1, 1, -1
       Ic = Ic * 2
       Ic = Ic + ie(j)
    end do

    if(Ic >= 2**nd) ier = 1

  end subroutine cld_ntrp_contract_indices


!...............................................................................
!...............................................................................
!> Cloud interpolation for linear cells
!> @param fvals values at the cell nodes
!> @param ts normalized [0,1]^nd cell coordinates
!> @param[out] interpolated value
!> @param[out] error flag, non zero if unsucessful
pure subroutine cld_ntrp_linear_cell_interp(fvals, ts, f, ier)
    real, intent(in) :: fvals(0:)
    real, intent(in) :: ts(:)
    real, intent(out):: f
    integer, intent(out) ::  ier

    integer j, nd, Ic, iflag
    integer ie(size(fvals))
    real    basis

    ier = 0
    f   = 0.
    nd   = size(ts)
    if(size(fvals) /= 2**nd) then
       ier = 1
       return
    endif

    do Ic = 0, 2**nd - 1
       basis = 1.
       call cld_ntrp_expand_index(Ic, ie, iflag)
       do j = 1, nd
          basis = basis * (  (1.0-real(ie(j)))*(1.0-ts(j)) + real(ie(j))*ts(j) )
       end do
       f = f + fvals(Ic)*basis
    end do

  end subroutine cld_ntrp_linear_cell_interp

!...............................................................................
!...............................................................................
pure subroutine cld_ntrp_locate_cell(axis, x, index, ier)
    real, intent(in)     :: axis(:) !< axis
    real, intent(in)     :: x       !< abscissae
    integer, intent(out) :: index   !< lower-left corner index
    integer, intent(out) ::  ier    !< error flag (0=ok)

    logical down
    integer n, index1, is
    real axis_1, axis_n, axis_min, axis_max
    ier   = 0
    index = -1
    down = .FALSE.
    n = size(axis)
    if(n < 2) then
       ier = 3
       return
    endif
    axis_1 = axis(1)
    axis_n = axis(n)
    axis_min = axis_1
    axis_max = axis_n
    if(axis_1 > axis_n) then
       down = .TRUE.
       axis_min = axis_n
       axis_max = axis_1
    endif

    if(x < axis_min-tol) then
       ier = 1
       return
    endif
    if(x > axis_max+tol) then
       ier = 2
       return
    endif

    index = floor(real(n-1)*(x - axis_1)/(axis_n-axis_1)) + 1
    index  = min(n-1, index)
    index1 = index+1

    if(.NOT. down) then
       if(axis(index) <= x+tol) then
          if(x <= axis(index1)+tol) then
             ! axis is uniform, or nearly so. Done!
             return
          else
             ! increase index
             is = index+1
             do index = is, n-1
                index1 = index+1
                if(axis(index1) >= x-tol) return
             enddo
          endif
       else
          ! decrease index
          is = index - 1
          do index = is, 1, -1
             if(axis(index) <= x+tol) return
          enddo
       endif
    else
       ! axis is pointing down
       if(axis(index) >= x-tol) then
          if(x >= axis(index1)-tol) then
             ! axis is uniform, or nearly so. Done!
             return
          else
             ! increase index
             is = index + 1
             do index = is, n-1
                index1 = index+1
                if(axis(index1) <= x+tol) return
             enddo
          endif
       else
          ! decrease index
          is = index - 1
          do index = is, 1, -1
             if(axis(index) >= x-tol) return
          enddo
       endif
    endif

  end subroutine cld_ntrp_locate_cell

!...............................................................................
!...............................................................................
pure subroutine cld_ntrp_get_flat_index(nsizes, indices, flat_index, ier)
    integer, intent(in)  :: nsizes(:)  !< size of array along each axis
    integer, intent(in)  :: indices(:) !< cell indices
    integer, intent(out) :: flat_index !< index into flattened array
    integer, intent(out) ::  ier       !< error flag (0=ok)

    integer nd, id

    ier = 0
    flat_index = -1
    nd = size(nsizes)
    if(nd /= size(indices)) then
       ! size mismatch
       ier = 1
       return
    endif

    flat_index = indices(nd)-1
    do id = nd-1, 1, -1
       flat_index = flat_index*nsizes(id) + indices(id)-1
    enddo
    flat_index = flat_index + 1

  end subroutine cld_ntrp_get_flat_index

!...............................................................................
!...............................................................................
pure subroutine cld_ntrp_get_cell_values(nsizes, fnodes, indices, fvals, ier)
    integer, intent(in)  :: nsizes(:)  !< size of fnodes along each axis
    real, intent(in)     :: fnodes(:)  !< flattened array of node values
    integer, intent(in)  :: indices(:) !< cell indices
    real, intent(out)    :: fvals(0:)  !< returned array values in the cell
    integer, intent(out) ::  ier       !< error flag (0=ok)

    integer id, nt, nd, flat_index, Ic, iflag
    integer, dimension(size(nsizes)) :: cell_indices, node_indices
    ier = 0
    fvals = 0.

    nd = size(nsizes)
    if(nd /= size(indices)) then
       ! size mismatch
       ier = 1
       return
    endif
    if(2**nd > size(fvals)) then
       ! not enough elements to hold result
       ier = 2
       return
    endif
    nt = 1
    do id = 1, nd
       nt = nt * nsizes(id)
    enddo
    if(nt /= size(fnodes)) then
       ! not enough node values
       ier = 3
       return
    endif

    do Ic = 0, 2**nd-1
       call cld_ntrp_expand_index(Ic, cell_indices, iflag)
       node_indices = indices + cell_indices
       call cld_ntrp_get_flat_index(nsizes, node_indices, flat_index, iflag)
       fvals(Ic) = fnodes(flat_index)
    enddo

  end subroutine cld_ntrp_get_cell_values

end MODULE cloud_interpolator_mod
!===============================================================================
!> @}
! close documentation grouping
