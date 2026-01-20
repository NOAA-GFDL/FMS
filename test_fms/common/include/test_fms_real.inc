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

  subroutine ARR_COMPARE_TOL_2D_ (arr0, arr1, tol, msg)
    real(FMS_TEST_KIND_), intent(in), dimension(:,:) :: arr0, arr1
    real(FMS_TEST_KIND_), intent(in) :: tol
    character(*), intent(in) :: msg

    if (any(abs(arr1 - arr0).gt.tol)) then
      call mpp_error(FATAL, "[2D] Unexpected result: " // msg)
    endif
  end subroutine ARR_COMPARE_TOL_2D_

  subroutine ARR_COMPARE_TOL_3D_ (arr0, arr1, tol, msg)
    real(FMS_TEST_KIND_), intent(in), dimension(:,:,:) :: arr0, arr1
    real(FMS_TEST_KIND_), intent(in) :: tol
    character(*), intent(in) :: msg

    if (any(abs(arr1 - arr0).gt.tol)) then
      call mpp_error(FATAL, "[3D] Unexpected result: " // msg)
    endif
  end subroutine ARR_COMPARE_TOL_3D_

  subroutine ARR_COMPARE_TOL_4D_ (arr0, arr1, tol, msg)
    real(FMS_TEST_KIND_), intent(in), dimension(:,:,:,:) :: arr0, arr1
    real(FMS_TEST_KIND_), intent(in) :: tol
    character(*), intent(in) :: msg

    if (any(abs(arr1 - arr0).gt.tol)) then
      call mpp_error(FATAL, "[4D] Unexpected result: " // msg)
    endif
  end subroutine ARR_COMPARE_TOL_4D_

  subroutine ARR_COMPARE_TOL_2D_SCALAR_ (arr, ans, tol, msg)
    real(FMS_TEST_KIND_), intent(in), dimension(:,:) :: arr
    real(FMS_TEST_KIND_), intent(in) :: ans
    real(FMS_TEST_KIND_), intent(in) :: tol
    character(*), intent(in) :: msg

    if (any(abs(arr - ans).gt.tol)) then
      call mpp_error(FATAL, "[2D] Unexpected result: " // msg)
    endif
  end subroutine ARR_COMPARE_TOL_2D_SCALAR_

  subroutine ARR_COMPARE_TOL_3D_SCALAR_ (arr, ans, tol, msg)
    real(FMS_TEST_KIND_), intent(in), dimension(:,:,:) :: arr
    real(FMS_TEST_KIND_), intent(in) :: ans
    real(FMS_TEST_KIND_), intent(in) :: tol
    character(*), intent(in) :: msg

    if (any(abs(arr - ans).gt.tol)) then
      call mpp_error(FATAL, "[3D] Unexpected result: " // msg)
    endif
  end subroutine ARR_COMPARE_TOL_3D_SCALAR_

  subroutine ARR_COMPARE_TOL_4D_SCALAR_ (arr, ans, tol, msg)
    real(FMS_TEST_KIND_), intent(in), dimension(:,:,:,:) :: arr
    real(FMS_TEST_KIND_), intent(in) :: ans
    real(FMS_TEST_KIND_), intent(in) :: tol
    character(*), intent(in) :: msg

    if (any(abs(arr - ans).gt.tol)) then
      call mpp_error(FATAL, "[4D] Unexpected result: " // msg)
    endif
  end subroutine ARR_COMPARE_TOL_4D_SCALAR_
