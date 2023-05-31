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

  !> Draws random scalar
  subroutine GET_RANDOM_NUMBER_0D_(stream, number)
    type(randomNumberStream), intent(inout) :: stream
    real(FMS_RN_KIND_),       intent(  out) :: number

    number = real(getRandomReal(stream%theNumbers), FMS_RN_KIND_)
  end subroutine

  !> Draws random 1D array
  subroutine GET_RANDOM_NUMBER_1D_(stream, numbers)
    type(randomNumberStream),         intent(inout) :: stream
    real(FMS_RN_KIND_), dimension(:), intent(  out) :: numbers

    ! Local variables
    integer :: i

    do i = 1, size(numbers)
      call GET_RANDOM_NUMBER_0D_(stream, numbers(i))
    end do
  end subroutine

  !> Draws random 2D array
  subroutine GET_RANDOM_NUMBER_2D_(stream, numbers)
    type(randomNumberStream),            intent(inout) :: stream
    real(FMS_RN_KIND_), dimension(:, :), intent(  out) :: numbers

    ! Local variables
    integer :: i

    do i = 1, size(numbers, 2)
      call GET_RANDOM_NUMBER_1D_(stream, numbers(:, i))
    end do
  end subroutine
