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

!> @brief Determines if a real input array has values which increase or
!!     decrease with strict monotonicity.
!! @return If the input array of real values either increases or decreases in a strictly monotonic manner,
!!     then true is returned. Otherwise, false is returned.

function MONOTONIC_ARRAY_(array, direction) result(ret)
real(FMS_MOD_KIND_), intent(in) :: array(:) !< An array of real values. If size(array) < 2, this function
                                            !! assumes the array is not monotonic; no fatal error will occur
                                            !! in this case.
integer, intent(out), optional :: direction !< If the input array is:
                                            !! >> strictly monotonic (small to large), then direction = +1.
                                            !! >> strictly monotonic (large to small), then direction = -1.
                                            !! >> not strictly monotonic, then direction = 0.
logical :: ret                              !< If the input array of real values either increases or
                                            !! decreases with strict monotonicity, then TRUE is returned;
                                            !! otherwise, FALSE is returned.
integer :: i

! initialize
  ret = .false.
  if (present(direction)) direction = 0

! array too short
  if ( size(array(:)) < 2 ) return

! ascending
  if ( array(1) < array(size(array(:))) ) then
     do i = 2, size(array(:))
       if (array(i-1) < array(i)) cycle
       return
     enddo
     ret = .true.
     if (present(direction)) direction = +1

! descending
  else
     do i = 2, size(array(:))
       if (array(i-1) > array(i)) cycle
       return
     enddo
     ret = .true.
     if (present(direction)) direction = -1
  endif

end function
