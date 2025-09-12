  subroutine ARR_INIT_2D_ (arr)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(out) :: arr(:,:)
    real(r8_kind) :: unif(size(arr,1), size(arr,2))
    type(randomNumberStream), save :: random_stream
    logical, save :: initialized = .false.

    if (.not.initialized) then
      random_stream = initializeRandomNumberStream(0)
      initialized = .true.
    endif

    call getRandomNumbers(random_stream, unif)

    arr = TYPECAST_ (1e9_r8_kind * (unif - 0.5_r8_kind), FMS_TEST_KIND_)
  end subroutine ARR_INIT_2D_

  subroutine ARR_INIT_3D_ (arr)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(out) :: arr(:,:,:)
    integer :: k

    do k = 1, size(arr, 3)
      call arr_init(arr(:, :, k))
    enddo
  end subroutine ARR_INIT_3D_

  subroutine ARR_COMPARE_2D_ (arr0, arr1, msg)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(in), dimension(:,:) :: arr0, arr1
    character(*), intent(in) :: msg

    if (any(arr0.ne.arr1)) then
      call mpp_error(FATAL, "Result from mpp_global_field (2D) does not agree with source data: " // msg)
    endif
  end subroutine ARR_COMPARE_2D_

  subroutine ARR_COMPARE_3D_ (arr0, arr1, msg)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(in), dimension(:,:,:) :: arr0, arr1
    character(*), intent(in) :: msg

    if (any(arr0.ne.arr1)) then
      call mpp_error(FATAL, "Result from mpp_global_field (3D) does not agree with source data: " // msg)
    endif
  end subroutine ARR_COMPARE_3D_
