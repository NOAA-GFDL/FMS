subroutine TEST_MONOTONIC_ARRAY_
  integer, parameter :: k = FMS_MOD_TEST_KIND_
  real(FMS_MOD_TEST_KIND_) :: arr1(1), arr2(2), arr5(5)

  arr1 = [1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false., 0)

  arr1 = [-1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false., 0)

  arr2 = [-1._k, 1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., 1)

  arr2 = [1._k, -1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., -1)

  arr2 = [1e10_k, 1e20_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., 1)

  arr2 = [1e10_k, 1e-20_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., -1)

  arr5 = [-2._k, -1._k, 0._k, 1._k, 2._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., 1)

  arr5 = [2._k, 1._k, 0._k, -1._k, -2._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., -1)

  arr5 = [1._k, 2._k, 3._k, 4._k, -5._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false., 0)

  arr5 = [-5._k, 4._k, 3._k, 2._k, 1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false., 0)

  arr5 = [1e-8_k, 1e-6_k, 1e-4_k, 1e-2_k, 1e0_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., 1)

  arr5 = [1e0_k, 1e-2_k, 1e-4_k, 1e-6_k, 1e-8_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., -1)

  arr5 = [1e10_k, 1e20_k, 1e30_k, 1e35_k, 9.99e37_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., 1)

  arr5 = [-1e10_k, -1e20_k, -1e30_k, -1e35_k, -9.99e37_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., -1)
end subroutine

subroutine TEST_MONOTONIC_ARRAY_ASSERT_(arr, ret_expected, direction_expected)
  real(FMS_MOD_TEST_KIND_), intent(in) :: arr(:)
  logical, intent(in) :: ret_expected
  integer, intent(in), optional :: direction_expected
  integer :: direction_test
  logical :: ret_test

  if (present(direction_expected)) then
    ret_test = monotonic_array(arr, direction_test)
    if (direction_test .ne. direction_expected) then
      write(stderr(), "(A)") "monotonic_array(arr, direction) returned incorrect direction: " // string(direction_test)
      write(stderr(), "(A)") "Expected value: " // string(direction_expected)
      call mpp_error(FATAL, "monotonic_array unit test failed")
    endif
  else
    ret_test = monotonic_array(arr)
  endif

  if (ret_test .neqv. ret_expected) then
    write(stderr(), "(A)") "monotonic_array(arr, direction) returned incorrect value: " // string(ret_test)
    write(stderr(), "(A)") "Expected value: " // string(ret_expected)
    call mpp_error(FATAL, "monotonic_array unit test failed")
  endif
end subroutine
