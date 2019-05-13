
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

      program test_quicksort
        implicit none
        integer :: list(16) = (/6, 2, 3, 4, 1, 45, 3432, 3245, 32545, 66555, 32, 1,3, -43254, 324, 54/)
        print *,'before list=', list
        call qksrt_quicksort(size(list), list, 1, size(list))
        print *,'after  list=', list
      end program test_quicksort
