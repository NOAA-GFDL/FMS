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

! Quicksort implementation
!   Simple implementation of a quicksort of an integer array.
!   Directly based on the pseudocode implementation at Wikipedia.
!   (https://en.wikipedia.org/w/index.php?title=Quicksort&oldid=883205571)
!   Notes:
!       - In-place sort; the original order is not preserved.
!       - Pivot is the middle value
!       - Hoare partitioning
!       - Unstable sort (like most Quicksort implementations)

module quicksort_mod
    implicit none
contains
    integer function hoare_partition(list, lo, hi) result(p)
        ! Hoare partition scheme
        !   - Iterate from each end towards middle, swap any inversions
        !   - Return index where iterations meet as partition point
        !   - Use median value as pivot test
        integer, intent(inout) :: list(:)
        integer, intent(in) :: lo, hi

        integer :: i, j, pivot, tmp
        integer :: i_new

        pivot = list((lo + hi) / 2)
        i = lo - 1
        j = hi + 1

        do
            i = i + 1
            do while (list(i) < pivot)
                i = i + 1
            end do

            j = j - 1
            do while (list(j) > pivot)
                j = j - 1
            end do

            if (i >= j) then
                exit
            end if

            tmp = list(i)
            list(i) = list(j)
            list(j) = tmp
        end do
        p = j
    end function hoare_partition

    recursive subroutine quicksort(list, lo_in, hi_in)
        integer, intent(inout) :: list(:)
        integer, intent(in), optional :: lo_in, hi_in

        integer :: lo, hi, p

        lo = 1; hi = size(list)
        if (present(lo_in)) lo = lo_in
        if (present(hi_in)) hi = hi_in

        if (lo < hi) then
            p = hoare_partition(list, lo, hi)
            call quicksort(list, lo, p)
            call quicksort(list, p + 1, hi)
        end if
    end subroutine quicksort
end module quicksort_mod
