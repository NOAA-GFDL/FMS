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
