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

!> @brief Utilities used in multiple test
module testing_utils
  use platform_mod,      only: r8_kind
  private

  public :: allocate_buffer
  public :: test_normal, test_openmp, test_halos
  public :: no_mask, logical_mask, real_mask

  integer, parameter :: test_normal = 0  !< sending a buffer in the compute domain
  integer, parameter :: test_openmp = 1  !< sending a buffer in the compute domain but with blocking
  integer, parameter :: test_halos = 2   !< sending a buffer in the data domain (i.e with halos)
  integer, parameter :: no_mask = 0      !< Not using a mask
  integer, parameter :: logical_mask = 1 !< Using a logical mask
  integer, parameter :: real_mask = 2    !< Using a real mask

  contains

  !> @brief Allocate the output buffer based on the starting/ending indices
  !! @return output buffer set to -999_r8_kind
  function allocate_buffer(is, ie, js, je, k, l) &
    result(buffer)
    integer, intent(in) :: is !< Starting x index
    integer, intent(in) :: ie !< Ending x index
    integer, intent(in) :: js !< Starting y index
    integer, intent(in) :: je !< Ending y index
    integer, intent(in) :: k  !< Number of points in the 4th dimension
    integer, intent(in) :: l  !< Number of points in the 5th dimension
    real(kind=r8_kind), allocatable :: buffer(:,:,:,:)

    allocate(buffer(is:ie, js:je, 1:k, 1:l))
    buffer = -999_r8_kind
  end function allocate_buffer
end module
