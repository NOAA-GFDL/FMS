#define C(x) x _ FMS_MOD_TEST_KIND

subroutine TEST_MONOTONIC_ARRAY
  real(FMS_MOD_TEST_KIND) :: arr1(1), arr2(2), arr5(5)

  arr1 = [C(1.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr1, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr1, .false., 0)

  arr1 = [C(-1.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr1, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr1, .false., 0)

  arr2 = [C(-1.), C(1.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true., 1)

  arr2 = [C(1.), C(-1.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true., -1)

  arr2 = [C(1e10), C(1e20)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true., 1)

  arr2 = [C(1e10), C(1e-20)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr2, .true., -1)

  arr5 = [C(-2.), C(-1.), C(0.), C(1.), C(2.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true., 1)

  arr5 = [C(2.), C(1.), C(0.), C(-1.), C(-2.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true., -1)

  arr5 = [C(1.), C(2.), C(3.), C(4.), C(-5.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .false., 0)

  arr5 = [C(-5.), C(4.), C(3.), C(2.), C(1.)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .false.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .false., 0)

  arr5 = [C(1e-8), C(1e-6), C(1e-4), C(1e-2), C(1e0)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true., 1)

  arr5 = [C(1e0), C(1e-2), C(1e-4), C(1e-6), C(1e-8)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true., -1)

  arr5 = [C(1e10), C(1e20), C(1e30), C(1e35), C(9.99e37)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true., 1)

  arr5 = [C(-1e10), C(-1e20), C(-1e30), C(-1e35), C(-9.99e37)]
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true.)
  call TEST_MONOTONIC_ARRAY_ASSERT(arr5, .true., -1)
end subroutine

subroutine TEST_MONOTONIC_ARRAY_ASSERT(arr, ret_expected, direction_expected)
  real(FMS_MOD_TEST_KIND), intent(in) :: arr(:)
  logical, intent(in) :: ret_expected
  integer, intent(in), optional :: direction_expected
  integer :: direction_test
  logical :: ret_test

  if (present(direction_expected)) then
    ret_test = monotonic_array(arr, direction_test)
    if (direction_test .ne. direction_expected) then
      write(stderr(), "(A)") "monotonic_array(arr, direction) returned incorrect direction: " // PRETTY(direction_test)
      write(stderr(), "(A)") "Expected value: " // PRETTY(direction_expected)
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

#undef C
