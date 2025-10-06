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
