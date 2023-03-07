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

!> @brief Determines if a real input array has monotonically increasing or
!!     decreasing values.
!! @return If the input array of real values either increases or decreases monotonically then true
!! is returned, otherwise false is returned.

function MONOTONIC_ARRAY_(array, direction)
real(FMS_MOD_KIND_), intent(in) :: array(:) !< An array of real values. If the size(array) < 2 this function
                                            !! assumes the array is not monotonic, no fatal error will occur.
integer, intent(out), optional :: direction !< If the input array is:
                                            !! >> monotonic (small to large) then direction = +1.
                                            !! >> monotonic (large to small) then direction = -1.
                                            !! >> not monotonic then direction = 0.
logical :: MONOTONIC_ARRAY_ !< If the input array of real values either increases or decreases monotonically
                            !! then TRUE is returned, otherwise FALSE is returned.
integer :: i

! initialize
  MONOTONIC_ARRAY_ = .false.
  if (present(direction)) direction = 0

! array too short
  if ( size(array(:)) < 2 ) return

! ascending
  if ( array(1) < array(size(array(:))) ) then
     do i = 2, size(array(:))
       if (array(i-1) < array(i)) cycle
       return
     enddo
     MONOTONIC_ARRAY_ = .true.
     if (present(direction)) direction = +1

! descending
  else
     do i = 2, size(array(:))
       if (array(i-1) > array(i)) cycle
       return
     enddo
     MONOTONIC_ARRAY_ = .true.
     if (present(direction)) direction = -1
  endif

end function
