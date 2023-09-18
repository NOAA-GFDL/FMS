subroutine TEST_MONOTONIC_ARRAY_
  integer, parameter :: k = FMS_MOD_TEST_KIND_
  real(FMS_MOD_TEST_KIND_) :: arr1(1), arr2(2), arr5(5)

  ! monotonic_array should return false when size=1

  arr1 = [1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false., 0)

  arr1 = [-1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr1, .false., 0)

  ! size=2, increasing
  arr2 = [-1._k, 1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., 1)

  ! size=2, decreasing
  arr2 = [1._k, -1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., -1)

  ! size=2, very large numbers, increasing
  arr2 = [1e10_k, 1e20_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., 1)

  ! size=2, very large numbers, decreasing
  arr2 = [1e10_k, 1e-20_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr2, .true., -1)

  ! Monotonically increasing, size=5
  arr5 = [-2._k, -1._k, 0._k, 1._k, 2._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., 1)

  ! Monotonically decreasing, size=5
  arr5 = [2._k, 1._k, 0._k, -1._k, -2._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., -1)

  ! Non-monotonic array, size=5
  arr5 = [1._k, 2._k, 3._k, 4._k, -5._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false., 0)

  ! Another permutation of non-monotonic array, size=5
  arr5 = [-5._k, 4._k, 3._k, 2._k, 1._k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false., 0)

  ! Monotonically increasing, size=5, small numbers
  arr5 = [1e-8_k, 1e-6_k, 1e-4_k, 1e-2_k, 1e0_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., 1)

  ! Monotonically decreasing, size=5, small numbers
  arr5 = [1e0_k, 1e-2_k, 1e-4_k, 1e-6_k, 1e-8_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., -1)

  ! Non-monotonic array, size=5, small numbers
  arr5 = [1e0_k, 1e-8_k, 1e-2_k, 1e-4_k, 1e-6_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false., 0)

  ! Monotonically increasing, size=5, positive large numbers
  arr5 = [1e10_k, 1e20_k, 1e30_k, 1e35_k, 9.99e37_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., 1)

  ! Monotonically decreasing, size=5, negative large numbers
  arr5 = [-1e10_k, -1e20_k, -1e30_k, -1e35_k, -9.99e37_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .true., -1)

  ! Non-monotonic array, size=5, negative large numbers
  arr5 = [-1e10_k, -1e20_k, -1e30_k, -9.99e37_k, -1e30_k]
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT_(arr5, .false., 0)
end subroutine

subroutine TEST_MONOTONIC_ARRAY_ASSERT_(arr, monotonic_expected, direction_expected)
  real(FMS_MOD_TEST_KIND_), intent(in) :: arr(:)
  logical, intent(in) :: monotonic_expected
  integer, intent(in), optional :: direction_expected
  integer :: direction_test
  logical :: monotonic_test

  if (present(direction_expected)) then
    monotonic_test = monotonic_array(arr, direction_test)
    if (direction_test .ne. direction_expected) then
      write(stderr(), "(A)") "monotonic_array(" // stringify(arr) // &
                           & ", direction) returned incorrect direction: " // string(direction_test)
      write(stderr(), "(A)") "Expected direction: " // string(direction_expected)
      call mpp_error(FATAL, "monotonic_array unit test failed")
    endif
  else
    monotonic_test = monotonic_array(arr)
  endif

  if (monotonic_test .neqv. monotonic_expected) then
    write(stderr(), "(A)") "monotonic_array(" // stringify(arr) // &
                         & ") returned incorrect value: " // string(monotonic_test)
    write(stderr(), "(A)") "Expected return value: " // string(monotonic_expected)
    call mpp_error(FATAL, "monotonic_array unit test failed")
  endif
end subroutine
