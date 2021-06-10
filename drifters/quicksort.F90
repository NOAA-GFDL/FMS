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
!> @cond
#undef _TYP
#define _TYP integer
!> @endcond

!> @file
!> @ingroup drifters
!> @brief Fortran implementation of quicksort
!> @author Magnus Lie Hetland

!> Create array partitions for quicksort
function qksrt_partition(n, list, start, end) result(top)
  implicit none
  integer, intent(in) :: n
  _TYP, intent(inout) :: list(n)
     integer, intent(in) :: start, end

     integer pivot, bottom, top
     logical done

     pivot = list(end)                          ! Partition around the last value
     bottom = start-1                           ! Start outside the area to be partitioned
     top = end                                  ! Ditto

     done = .false.
     do while (.not. done)                      ! Until all elements are partitioned...

        do while (.not. done)                  ! Until we find an out of place element...
           bottom = bottom+1                  ! ... move the bottom up.

           if(bottom == top) then             ! If we hit the top...
              done = .true.                  ! ... we are done.
              exit
           endif

           if(list(bottom) > pivot) then           ! Is the bottom out of place?
              list(top) = list(bottom)       ! Then put it at the top...
              exit                          ! ... and start searching from the top.
           endif
        enddo

        do while (.not. done)                        ! Until we find an out of place element...
           top = top-1                        ! ... move the top down.

           if(top == bottom) then                  ! If we hit the bottom...
              done = .true.                      ! ... we are done.
              exit
           endif

           if(list(top) < pivot) then              ! Is the top out of place?
              list(bottom) = list(top)       ! Then put it at the bottom...
              exit                          ! ...and start searching from the bottom.
           endif
        enddo
     enddo

     list(top) = pivot                          ! Put the pivot in its place.
     ! Return the split point

end function qksrt_partition

!> quicksort a given list
recursive subroutine qksrt_quicksort(n, list, start, end)
     implicit none
     integer, intent(in) :: n
     _TYP, intent(inout) :: list(n)
     integer, intent(in) :: start, end
     integer :: split, qksrt_partition
     external :: qksrt_partition
     if(start < end) then                            ! If there are two or more elements...
        split = qksrt_partition(n, list, start, end)    ! ... partition the sublist...
        call qksrt_quicksort(n, list,  start, split-1)        ! ... and sort both halves.
        call qksrt_quicksort(n, list, split+1, end)
     endif
end subroutine qksrt_quicksort
