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

program test_cloud_interpolator
  use cloud_interpolator_mod
  use mpp_mod, only : mpp_error, FATAL, stdout, mpp_init, mpp_exit

  implicit none

  call mpp_init()

  call test_expansion_contraction
  call test_linear_cell_interpolation
  call test_cell_search
  call test_get_node_values

  call mpp_exit()
CONTAINS
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
         call mpp_error(FATAL,'ERROR: contraction/expansion test failed (ie1/=ie2)')
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
      if(ier/=0) then
         print *,'ERROR flag ier=', ier
         call mpp_error(FATAL,'ERROR: linear cell interpolation test failed (ier/=0)')
      endif
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

      if(ier_tot /= 0) then
         print *,'ERROR flag ier_tot=', ier_tot
         call mpp_error(FATAL,'ERROR: cell search test failed (ier_tot/=0)')
      endif

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
      call cld_ntrp_get_cell_values((/n1,n2,n3/), pack(fnodes, .TRUE.) , indices, fvals, ier)
      fexact = (/0.0, 0.2, 0.0, 0.2, 0.0, 0.2, 0.5, 0.7/)
      if(ier/=0) print *,'ERROR flag ier=', ier
      print *,'indices ', indices
      print *,'fvals=', fvals, ' ==? ', fexact
      error_tot = error_tot + abs(sum(fvals - fexact))

      indices = (/5,4,2/)
      call cld_ntrp_get_cell_values((/n1,n2,n3/), pack(fnodes, .TRUE.) , indices, fvals, ier)
      fexact = (/2.3, 2.5, 2.8, 3.0, 6.8, 7.0, 8.8, 9.0/)
      if(ier/=0) print *,'ERROR flag ier=', ier
      print *,'indices ', indices
      print *,'fvals=', fvals, ' ==? ', fexact
      error_tot = error_tot + abs(sum(fvals - fexact))

      if(error_tot /= 0) then
         print *,'ERROR flag error_tot=', error_tot
         call mpp_error(FATAL,'ERROR: get node values test failed (error_tot/=0)')
      endif

end subroutine test_get_node_values

  end program test_cloud_interpolator
