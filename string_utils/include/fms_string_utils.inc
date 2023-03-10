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

!> @brief Converts a 1D array of real numbers to a string
!> @return The 1D array as a string
function STRINGIFY_1D_(arr, fmt)
  real(STRING_UTILS_KIND_), dimension(:), intent(in) :: arr !< Real array to be converted to a string
  character(*), intent(in), optional :: fmt !< Optional format string for the real array entries
  character(:), allocatable :: STRINGIFY_1D_
  integer :: i, n

  n = size(arr)

  if (n .gt. 0) then
    STRINGIFY_1D_ = "[" // string(arr(1), fmt)
  else
    STRINGIFY_1D_ = "["
  endif

  do i = 2,n
    STRINGIFY_1D_ = STRINGIFY_1D_ // ", " // string(arr(i), fmt)
  enddo

  STRINGIFY_1D_ = STRINGIFY_1D_ // "]"
end function

!> @brief Converts a 2D array of real numbers to a string
!> @return The 2D array as a string
function STRINGIFY_2D_(arr, fmt)
  real(STRING_UTILS_KIND_), dimension(:,:), intent(in) :: arr !< Real array to be converted to a string
  character(*), intent(in), optional :: fmt !< Optional format string for the real array entries
  character(:), allocatable :: STRINGIFY_2D_
  integer :: i, n

  n = size(arr, 2)

  if (n .gt. 0) then
    STRINGIFY_2D_ = "[" // STRINGIFY_1D_(arr(:,1), fmt)
  else
    STRINGIFY_2D_ = "["
  endif

  do i = 2,n
    STRINGIFY_2D_ = STRINGIFY_2D_ // ", " // STRINGIFY_1D_(arr(:,i), fmt)
  enddo

  STRINGIFY_2D_ = STRINGIFY_2D_ // "]"
end function

!> @brief Converts a 3D array of real numbers to a string
!> @return The 3D array as a string
function STRINGIFY_3D_(arr, fmt)
  real(STRING_UTILS_KIND_), dimension(:,:,:), intent(in) :: arr !< Real array to be converted to a string
  character(*), intent(in), optional :: fmt !< Optional format string for the real array entries
  character(:), allocatable :: STRINGIFY_3D_
  integer :: i, n

  n = size(arr, 3)

  if (n .gt. 0) then
    STRINGIFY_3D_ = "[" // STRINGIFY_2D_(arr(:,:,1), fmt)
  else
    STRINGIFY_3D_ = "["
  endif

  do i = 2,n
    STRINGIFY_3D_ = STRINGIFY_3D_ // ", " // STRINGIFY_2D_(arr(:,:,i), fmt)
  enddo

  STRINGIFY_3D_ = STRINGIFY_3D_ // "]"
end function
