!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
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
