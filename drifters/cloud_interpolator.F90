! nf95 -r8 -g -I ~/regression/ia64/23-Jun-2005/CM2.1U_Control-1990_E1.k32pe/include/ -D_TEST_CLOUD_INTERPOLATOR -D_F95 cloud_interpolator.F90

#include <fms_platform.h>

#define _FLATTEN(A) reshape((A), (/size((A))/) )

MODULE cloud_interpolator_mod
  implicit none
  private

  public :: cld_ntrp_linear_cell_interp, cld_ntrp_locate_cell, cld_ntrp_get_cell_values
#ifdef _TEST_CLOUD_INTERPOLATOR
  public :: cld_ntrp_expand_index, cld_ntrp_contract_indices
#endif

! Include variable "version" to be written to log file.
#include<file_version.h>
real, parameter           :: tol = 10.0*epsilon(1.)

CONTAINS

!...............................................................................
  _PURE subroutine cld_ntrp_expand_index(Ic, ie, ier)
    integer, intent(in)  ::  Ic    ! contacted index
    integer, intent(out) ::  ie(:) ! expanded list of indices
    integer, intent(out) ::  ier   ! error flag (0=ok)

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
!...............................................................................
  _PURE subroutine cld_ntrp_contract_indices(ie, Ic, ier)
    integer, intent(in) ::  ie(:)  ! expanded list of indices
    integer, intent(out)  ::  Ic   ! contacted index
    integer, intent(out) ::  ier   ! error flag (0=ok)

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
  _PURE subroutine cld_ntrp_linear_cell_interp(fvals, ts, f, ier)
    real, intent(in) :: fvals(0:)  ! values at the cell nodes
    real, intent(in) :: ts(:)      ! normalized [0,1]^nd cell coordinates
    real, intent(out):: f          ! interpolated value
    integer, intent(out) ::  ier   ! error flag (0=ok)
    
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
  _PURE subroutine cld_ntrp_locate_cell(axis, x, index, ier)
    real, intent(in)     :: axis(:) ! axis 
    real, intent(in)     :: x       ! abscissae
    integer, intent(out) :: index   ! lower-left corner index
    integer, intent(out) ::  ier    ! error flag (0=ok)

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
  _PURE subroutine cld_ntrp_get_flat_index(nsizes, indices, flat_index, ier)
    integer, intent(in)  :: nsizes(:)  ! size of array along each axis
    integer, intent(in)  :: indices(:) ! cell indices
    integer, intent(out) :: flat_index ! index into flattened array
    integer, intent(out) ::  ier       ! error flag (0=ok)

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
  _PURE subroutine cld_ntrp_get_cell_values(nsizes, fnodes, indices, fvals, ier)
    integer, intent(in)  :: nsizes(:)  ! size of fnodes along each axis
    real, intent(in)     :: fnodes(:)  ! flattened array of node values
    integer, intent(in)  :: indices(:) ! cell indices
    real, intent(out)    :: fvals(0:)  ! returned array values in the cell
    integer, intent(out) ::  ier       ! error flag (0=ok)

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

#ifdef _TEST_CLOUD_INTERPOLATOR
program test
  use cloud_interpolator_mod
  implicit none

  call test_expansion_contraction
  call test_linear_cell_interpolation
  call test_cell_search
  call test_get_node_values

  contains
    subroutine test_expansion_contraction
      integer ie1(4), ie2(4), Ic, ier, idiff, j
      ie1 = (/1,0,1,1/)
      call cld_ntrp_contract_indices(ie1, Ic, ier)
      if(ier/=0) print *,'ERROR flag ier=', ier
      call cld_ntrp_expand_index(Ic, ie2, ier)
      if(ier/=0) print *,'ERROR flag ier=', ier
      idiff = 0
      do j = 1, size(ie1)
         idiff = idiff + abs(ie1(j)-ie2(j))
      end do
      if(idiff/=0) then
         print *,'ERROR: contraction/expansion test failed (ie1/=ie2)'
      endif
      print *,'ie1 = ', ie1
      print *,'ie2 = ', ie2
      print *,'Ic  = ', Ic
      
    end subroutine test_expansion_contraction

    subroutine test_linear_cell_interpolation
      integer, parameter :: nd = 3 
      real :: fvals(2**nd), ts(nd)
      real :: fi, fx
      integer ier
      ! f  = 1 + x + 2*y + 3*z
      fvals = (/ 1., 2., 3., 4., 4., 5., 6., 7. /)
      ts    = (/ 0.1, 0.2, 0.3 /)
      fx    = 1. + ts(1) + 2*ts(2) + 3*ts(3)
      call cld_ntrp_linear_cell_interp(fvals, ts, fi, ier)
      if(ier/=0) print *,'ERROR flag ier=', ier
      print *,'fi, fx = ', fi, fx
    end subroutine test_linear_cell_interpolation

    subroutine test_cell_search
      integer index, ier
      integer, parameter :: n = 5
      real :: axis1(n) = (/0., 0.1, 0.2, 0.3, 0.4/)
      real :: axis2(n) = (/0., 0.01, 0.02, 0.03, 0.4/)
      real :: axis3(n) = (/0.4, 0.3, 0.2, 0.1, 0./)
      real :: axis4(n) = (/0.4, 0.03, 0.02, 0.01, 0./)
      real x
      integer :: ier_tot = 0
      
      print *,'axis1=', axis1
      x = -0.0001
      call cld_ntrp_locate_cell(axis1, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))
      x = 0.
      call cld_ntrp_locate_cell(axis1, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 1
      ier_tot = ier_tot + abs(index - (1))
      x = 0.1
      call cld_ntrp_locate_cell(axis1, x, index, ier)
      print *, ' x=',x, ' index=', index, ' ==? ', 2
      ier_tot = ier_tot + abs(index - (2))
      x = 0.4
      call cld_ntrp_locate_cell(axis1, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 4
      ier_tot = ier_tot + abs(index - (4))
      x = 0.40001
      call cld_ntrp_locate_cell(axis1, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))

      print *,'axis2=', axis1
      x = -0.0001
      call cld_ntrp_locate_cell(axis2, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))
      x = 0.
      call cld_ntrp_locate_cell(axis2, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 1
      ier_tot = ier_tot + abs(index - (1))
      x = 0.1
      call cld_ntrp_locate_cell(axis2, x, index, ier)
      print *, ' x=',x, ' index=', index, ' ==? ', 4
      ier_tot = ier_tot + abs(index - (4))
      x = 0.4
      call cld_ntrp_locate_cell(axis2, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 4
      ier_tot = ier_tot + abs(index - (4))
      x = 0.40001
      call cld_ntrp_locate_cell(axis2, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))

      print *,'axis3=', axis1
      x = -0.0001
      call cld_ntrp_locate_cell(axis3, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))
      x = 0.
      call cld_ntrp_locate_cell(axis3, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 4
      ier_tot = ier_tot + abs(index - (4))
      x = 0.1
      call cld_ntrp_locate_cell(axis3, x, index, ier)
      print *, ' x=',x, ' index=', index, ' ==? ', 4
      ier_tot = ier_tot + abs(index - (4))
      x = 0.4
      call cld_ntrp_locate_cell(axis3, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 1
      ier_tot = ier_tot + abs(index - (1))
      x = 0.40001
      call cld_ntrp_locate_cell(axis3, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))

      print *,'axis4=', axis1
      x = -0.0001
      call cld_ntrp_locate_cell(axis4, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))
      x = 0.
      call cld_ntrp_locate_cell(axis4, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 4
      ier_tot = ier_tot + abs(index - (4))
      x = 0.1
      call cld_ntrp_locate_cell(axis4, x, index, ier)
      print *, ' x=',x, ' index=', index, ' ==? ', 1
      ier_tot = ier_tot + abs(index - (1))
      x = 0.4
      call cld_ntrp_locate_cell(axis4, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', 1
      ier_tot = ier_tot + abs(index - (1))
      x = 0.40001
      call cld_ntrp_locate_cell(axis4, x, index, ier)
      print *,' x=',x, ' index=', index, ' ==? ', -1
      ier_tot = ier_tot + abs(index - (-1))

      print *,'Total error in test_cell_search: ', ier_tot

    end subroutine test_cell_search

    subroutine test_get_node_values
      integer, parameter :: nd = 3, n1=6, n2=5, n3=4
      real, dimension(n1, n2, n3) :: fnodes
      real :: fvals(2**nd), fexact(2**nd)
      real x, y, z
      integer i, j, k, ier, indices(nd)
      real :: error_tot = 0.
      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               x = 1* real(i-1)/real(n1-1)
               y = 2* real(j-1)/real(n2-1)
               z = 3* real(k-1)/real(n3-1)
               fnodes(i,j,k) = x + y*z**2
            enddo
         enddo
      enddo

      indices = (/1,1,1/)
      call cld_ntrp_get_cell_values((/n1,n2,n3/), _FLATTEN(fnodes), indices, fvals, ier)
      fexact = (/0.0, 0.2, 0.0, 0.2, 0.0, 0.2, 0.5, 0.7/)
      if(ier/=0) print *,'ERROR flag ier=', ier
      print *,'indices ', indices
      print *,'fvals=', fvals, ' ==? ', fexact
      error_tot = error_tot + abs(sum(fvals - fexact))

      indices = (/5,4,2/)
      call cld_ntrp_get_cell_values((/n1,n2,n3/), _FLATTEN(fnodes), indices, fvals, ier)
      fexact = (/2.3, 2.5, 2.8, 3.0, 6.8, 7.0, 8.8, 9.0/)
      if(ier/=0) print *,'ERROR flag ier=', ier
      print *,'indices ', indices
      print *,'fvals=', fvals, ' ==? ', fexact
      error_tot = error_tot + abs(sum(fvals - fexact))

      print *,'Total error in test_get_node_values: ', error_tot


    end subroutine test_get_node_values

  end program test

#endif
