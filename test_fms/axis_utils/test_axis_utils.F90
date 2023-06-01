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

program test_axis_utils

use fms_mod,         only : fms_init, fms_end, lowercase
use fms2_io_mod, only: FmsNetcdfFile_t, open_file, close_file, register_axis, register_field, &
                     & register_variable_attribute, write_data
use platform_mod, only: r4_kind, r8_kind
use mpp_mod, only: mpp_error, fatal, stderr
use fms_string_utils_mod, only: string, stringify
use axis_utils2_mod

implicit none

type GetAxisCartTest_t
  type(FmsNetcdfFile_t) :: fileobj
  type(GetAxisCartTestCase_t), pointer :: test0, test1
end type

type GetAxisCartTestCase_t
  character(:), allocatable :: var
  character(1) :: cart
  type(GetAxisCartTestCase_t), pointer :: next => NULL()
end type

integer, parameter :: k = AU_TEST_KIND_
real(k), parameter :: pi = 4._k * atan(1._k)

integer :: i
character(100) :: arg

call fms_init

do i=1,command_argument_count()
  call get_command_argument(i, arg)

  select case (arg)
    case ('--get-axis-modulo')
      print "(A)", "Testing get_axis_modulo"
      call test_get_axis_modulo

    case ('--get-axis-modulo-times')
      print "(A)", "Testing get_axis_modulo_times"
      call test_get_axis_modulo_times

    case ('--get-axis-cart')
      print "(A)", "Testing get_axis_cart"
      call test_get_axis_cart

    case ('--lon-in-range')
      print "(A)", "Testing lon_in_range"
      call test_lon_in_range

    case ('--frac-index')
      print "(A)", "Testing frac_index"
      call test_frac_index

    case ('--frac-index-fail')
      print "(A)", "Testing frac_index (FAILURE)"
      call test_frac_index_fail

    case ('--nearest-index-increasing')
      print "(A)", "Testing nearest_index with a monotonically increasing array"
      call test_nearest_index(.true.)

    case ('--nearest-index-decreasing')
      print "(A)", "Testing nearest_index with a monotonically decreasing array"
      call test_nearest_index(.false.)

    case ('--nearest-index-fail')
      print "(A)", "Testing nearest_index (FAILURE)"
      call test_nearest_index_fail

    case ('--axis-edges')
      print "(A)", "Testing axis_edges"
      call test_axis_edges

    case ('--tranlon')
      print "(A)", "Testing tranlon"
      call test_tranlon

    case ('--interp-1d-1d')
      print "(A)", "Testing interp_1d_1d"
      call test_interp_1d_1d

    case ('--interp-1d-2d')
      print "(A)", "Testing interp_1d_2d"
      call test_interp_1d_2d

    case ('--interp-1d-3d')
      print "(A)", "Testing interp_1d_3d"
      call test_interp_1d_3d

    case default
      write(stderr(),"(A)") "Unrecognized command line option: " // trim(arg)
  end select
enddo

call fms_end

contains

! Status: TODO
! function get_axis_modulo(fileobj, axisname)
subroutine test_get_axis_modulo
  type(FmsNetcdfFile_t) :: fileobj

  write(stderr(), "(A)") "Warning: get_axis_modulo unit test not yet implemented"
end subroutine

! Status: TODO
! function get_axis_modulo_times(fileobj, axisname, tbeg, tend)
subroutine test_get_axis_modulo_times
  type(FmsNetcdfFile_t) :: fileobj

  write(stderr(), "(A)") "Warning: get_axis_modulo_times unit test not yet implemented"
end subroutine

subroutine test_get_axis_cart
  type(GetAxisCartTest_t) :: test
  type(GetAxisCartTestCase_t), pointer :: test_nonexistent_var
  character(:), allocatable :: var_name, attr_name, attr_value
  integer :: i, j

  character(*), parameter, dimension(*) :: &
    & special_axis_names_x = [character(12) :: "lon", "x", "degrees_e", "degrees_east", "degreese"], &
    & special_axis_names_y = [character(13) :: "lat", "y", "degrees_n", "degrees_north", "degreesn"], &
    & special_axis_names_z = [character(6) :: "depth", "height", "z", "cm", "m", "pa", "hpa"], &
    & special_axis_names_t = [character(4) :: "time", "t", "sec", "min", "hou", "day", "mon", "yea"], &
    & attr_names           = [character(14) :: "cartesian_axis", "axis"], &
    & xyzt_uc              = ["X", "Y", "Z", "T"]

  call open_netcdf_w(test%fileobj)
  call register_axis(test%fileobj, "dim1", 1)

  ! Check a variable which does not exist

  allocate(test_nonexistent_var)
  test_nonexistent_var%var = "does_not_exist"
  test_nonexistent_var%cart = "N"

  test%test0 => test_nonexistent_var
  test%test1 => test_nonexistent_var

  ! Check a variable which exists, but which has neither a "cartesian_axis" nor an "axis" attribute.
  var_name = "exists_no_attributes"
  call get_axis_cart_test_add(test, var_name, "N")

  do i=1,size(attr_names)
    attr_name = trim(attr_names(i))

    ! Check an unknown value on a "cartesian_axis" or "axis" attribute.
    ! TODO: This test fails. It should be uncommented if/when get_axis_cart's behavior is fixed.

    !attr_value = "unexpected"
    !var_name = attr_name // "_attr_value_" // attr_value
    !call get_axis_cart_test_add(test, var_name, "N")
    !call register_variable_attribute(test%fileobj, var_name, attr_name, attr_value, str_len=len(attr_value))

    do j=1,size(xyzt_uc)
      ! Check upper-case "axis" attributes"
      attr_value = xyzt_uc(j)
      var_name = attr_name // "_attr_value_" // attr_value
      call get_axis_cart_test_add(test, var_name, xyzt_uc(j))
      call register_variable_attribute(test%fileobj, var_name, attr_name, attr_value, str_len=len(attr_value))

      ! Check lower-case "axis" attributes"
      attr_value = lowercase(xyzt_uc(j))
      var_name = attr_name // "_attr_value_" // attr_value
      call get_axis_cart_test_add(test, var_name, xyzt_uc(j))
      call register_variable_attribute(test%fileobj, var_name, attr_name, attr_value, str_len=len(attr_value))
    enddo
  enddo

  call test_special_axis_names(test, special_axis_names_x, "X")
  call test_special_axis_names(test, special_axis_names_y, "Y")
  call test_special_axis_names(test, special_axis_names_z, "Z")
  call test_special_axis_names(test, special_axis_names_t, "T")

  call close_file(test%fileobj)

  call get_axis_cart_tests_run(test)
end subroutine

subroutine get_axis_cart_test_add(test, var_name, cart)
  type(GetAxisCartTest_t), intent(inout) :: test
  type(GetAxisCartTestCase_t), pointer :: test_case
  character(*), intent(in) :: var_name
  character(1), intent(in) :: cart
  character(:), allocatable :: kind_str

  if (k .eq. r4_kind) then
    kind_str = "float"
  else
    kind_str = "double"
  endif

  call register_field(test%fileobj, var_name, kind_str, dimensions=["dim1"])

  allocate(test_case)
  test_case%var = var_name
  test_case%cart = cart

  test%test1%next => test_case
  test%test1 => test_case
end subroutine

subroutine get_axis_cart_tests_run(test)
  type(GetAxisCartTest_t), intent(inout) :: test
  type(GetAxisCartTestCase_t), pointer :: test_case, next
  character(1) :: cart_test
  integer :: i

  call open_netcdf_r(test%fileobj)

  test_case => test%test0

  do while (associated(test_case))
    cart_test = " "
    call get_axis_cart(test%fileobj, test_case%var, cart_test)

    if (cart_test .ne. test_case%cart) then
      write(stderr(), "(A)") "get_axis_cart result for variable '" // test_case%var // "': " // cart_test
      write(stderr(), "(A)") "Expected result: " // test_case%cart
      call mpp_error(FATAL, "get_axis_cart unit test failed")
    endif

    next => test_case%next
    deallocate(test_case)
    test_case => next
  enddo

  call close_file(test%fileobj)
end subroutine

subroutine test_special_axis_names(test, special_axis_names, ret_expected)
  type(GetAxisCartTest_t), intent(inout) :: test
  character(*), intent(in) :: special_axis_names(:), ret_expected
  character(:), allocatable :: var_name
  integer :: i

  do i=1,size(special_axis_names)
    var_name = trim(special_axis_names(i))
    call get_axis_cart_test_add(test, var_name, ret_expected)
  enddo
end subroutine

subroutine test_lon_in_range
  real(k), parameter :: eps_big = 1e-3_k, eps_tiny = 1e-5_k
  real(k), parameter :: pi_plus_360 = 360._k + pi

  ! Test some cases where no translation is needed
  call lon_in_range_assert(0._k,              0._k,  0._k)
  call lon_in_range_assert(1._k,              0._k,  1._k)
  call lon_in_range_assert(350._k,            0._k,  350._k)
  call lon_in_range_assert(1._k,              1._k,  1._k)
  call lon_in_range_assert(350._k,            1._k,  350._k)
  call lon_in_range_assert(359._k,            0._k,  359._k)
  call lon_in_range_assert(359._k,            1._k,  359._k)
  call lon_in_range_assert(pi,                0._k,  pi)

  ! Test up-translation
  call lon_in_range_assert(-2._k,             -1._k, 358._k)
  call lon_in_range_assert(-2._k,             0._k,  358._k)
  call lon_in_range_assert(-2._k,             5._k,  358._k)
  call lon_in_range_assert(-1._k,             0._k,  359._k)
  call lon_in_range_assert(-1._k,             5._k,  359._k)
  call lon_in_range_assert(0._k,              5._k,  360._k)
  call lon_in_range_assert(1._k,              5._k,  361._k)
  call lon_in_range_assert(-pi,               0._k,  360._k - pi)

  ! Test down-translation
  call lon_in_range_assert(359._k,            -1._k, -1._k)
  call lon_in_range_assert(360._k,            -1._k, 0._k)
  call lon_in_range_assert(360._k,            0._k,  0._k)
  call lon_in_range_assert(361._k,            -1._k, 1._k)
  call lon_in_range_assert(361._k,            0._k,  1._k)
  call lon_in_range_assert(362._k,            -1._k, 2._k)
  call lon_in_range_assert(362._k,            0._k,  2._k)
  call lon_in_range_assert(pi_plus_360,       0._k,  pi_plus_360 - 360._k)

  ! Test rounding behavior
  call lon_in_range_assert(eps_tiny,          0._k,  0._k)
  call lon_in_range_assert(eps_big,           0._k,  eps_big)
  call lon_in_range_assert(360._k - eps_tiny, 0._k,  0._k)
  call lon_in_range_assert(360._k - eps_big,  0._k,  360._k - eps_big)
end subroutine

subroutine lon_in_range_assert(lon, l_start, ret_expected)
  real(k), intent(in) :: lon, l_start, ret_expected
  real(k) :: ret_test

  ret_test = lon_in_range(lon, l_start)

  if (ret_test /= ret_expected) then
    write(stderr(), "(A)") "lon_in_range(" // string(lon) // ", " // string(l_start) // &
                         & ") returned erroneous value: " // string(ret_test)
    write(stderr(), "(A)") "Expected return value: " // string(ret_expected)
    call mpp_error(FATAL, "lon_in_range unit test failed")
  endif
end subroutine

#define CALC_FRAC_INDEX_(i, v, values) real(i, k) + (v - values(i)) / (values(i + 1) - values(i))

subroutine test_frac_index
  real(k) :: values(6), v, fi
  integer :: i, n
  real(k), parameter :: f10=.1_k, f25=.25_k, f50=.5_k, f99=.99_k

  values = [1._k, 2._k, 3._k, 5._k, 10._k, 11._k]
  n = size(values)

  ! Test values outside of the input array
  call frac_index_assert(real(values(1), k) - f50, values, -1._k)
  call frac_index_assert(real(values(n), k) + f50, values, -1._k)

  ! Test the actual indices
  do i=1,n
    v = values(i)
    call frac_index_assert(v, values, real(i, k))
  enddo

  ! Test the 10% point
  do i=1,n-1
    v = values(i) + f10*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX_(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo

  ! Test the 25% point
  do i=1,n-1
    v = values(i) + f25*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX_(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo

  ! Test the mid-point
  do i=1,n-1
    v = values(i) + f50*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX_(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo

  ! Test the 99% point
  do i=1,n-1
    v = values(i) + f99*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX_(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo
end subroutine

subroutine frac_index_assert(fval, arr, ret_expected)
  real(k), intent(in) :: fval, arr(:), ret_expected
  real(k) :: ret_test

  ret_test = frac_index(fval, arr)

  if (ret_test /= ret_expected) then
    write(stderr(), "(A)") "frac_index(" // string(fval) // ", " // stringify(arr) // &
                         & ") returned erroneous value: " // string(ret_test)
    write(stderr(), "(A)") "Expected return value: " // string(ret_expected)
    call mpp_error(FATAL, "frac_index unit test failed")
  endif
end subroutine

! Test that frac_index fails with a non-monotonic array
subroutine test_frac_index_fail
  real(k) :: values(5)
  real(k) :: ret_test

  values = [1._k, 2._k, 4._k, 3._k, 5._k]
  ret_test = frac_index(1.5_k, values)
end subroutine

subroutine test_nearest_index(increasing_array)
  logical, intent(in) :: increasing_array !< .True. if test using an increasing array
  real(k) :: arr(5)
  integer :: ans(12)

  arr = [5._k, 12._k, 20._k, 40._k, 100._k]
  if (increasing_array) then
    ans=(/1, 5, 1, 2, 3, 4, 5, 1, 2, 2, 3, 3/)
  else
    arr = arr(ubound(arr,dim=1)::-1)
    ans=(/5, 1, 5, 4, 3, 2, 1, 5, 4, 4, 3, 3/)
  endif

  ! Test values beyond array boundaries
  call nearest_index_assert(4._k,    arr, ans(1))
  call nearest_index_assert(1000._k, arr, ans(2))

  ! Test values actually in the array
  call nearest_index_assert(5._k,    arr, ans(3))
  call nearest_index_assert(12._k,   arr, ans(4))
  call nearest_index_assert(20._k,   arr, ans(5))
  call nearest_index_assert(40._k,   arr, ans(6))
  call nearest_index_assert(100._k,  arr, ans(7))

  ! Test the intervals between array values
  call nearest_index_assert(6._k,    arr, ans(8))
  call nearest_index_assert(11._k,   arr, ans(9))
  call nearest_index_assert(15._k,   arr, ans(10))
  call nearest_index_assert(18._k,   arr, ans(11))
  call nearest_index_assert(29._k,   arr, ans(12))
end subroutine

subroutine nearest_index_assert(val, arr, ret_expected)
  real(k), intent(in) :: val, arr(:)
  integer, intent(in) :: ret_expected
  integer :: ret_test

  ret_test = nearest_index(val, arr)

  if (ret_test /= ret_expected) then
    write(stderr(), "(A)") "nearest_index(" // string(val) // ", " // stringify(arr) // &
                         & ") returned erroneous value: " // string(ret_test)
    write(stderr(), "(A)") "Expected return value: " // string(ret_expected)
    call mpp_error(FATAL, "nearest_index unit test failed")
  endif
end subroutine

! Test that nearest_index fails with a non-monotonic array
subroutine test_nearest_index_fail
  real(k) :: arr(5)
  integer :: ret_test

  arr=[5._k, 12._k, 40._k, 20._k, 100._k]
  ret_test = nearest_index(5._k, arr)
end subroutine

subroutine test_axis_edges
  real(k) :: data_in_var(10)
  real(k) :: data_in_var_edges(2,10)
  real(k) :: data_in_answers(11)
  type(FmsNetcdfFile_t) :: fileobj
  real(k)    :: answers(11)
  integer :: i

  do i=1,10
     data_in_var(i) = real(i, k) - 0.5_k

     data_in_var_edges(1,i) = real(i-1, k)
     data_in_var_edges(2,i) = real(i, k)

     data_in_answers(i) = real(i-1, k)
  enddo

  data_in_answers(11) = 10._k

  call open_netcdf_w(fileobj)

  call register_axis(fileobj, "dim1", 10)
  call register_axis(fileobj, "dim2", 2)

  call register_field(fileobj, "axis", "double", dimensions=["dim1"])

  call register_field(fileobj, "axis_with_bounds", "double", dimensions=["dim1"])
  call register_variable_attribute(fileobj, "axis_with_bounds", "bounds", "bounds", str_len=6)
  call register_field(fileobj, "bounds", "double", dimensions=["dim2", "dim1"])

  call register_field(fileobj, "axis_with_edges", "double", dimensions=["dim1"])
  call register_variable_attribute(fileobj, "axis_with_edges", "edges", "edges"//char(0), str_len=6)
  call register_field(fileobj, "edges", "double", dimensions=["dim2", "dim1"])

  call write_data(fileobj, "axis", data_in_var)
  call write_data(fileobj, "axis_with_bounds", data_in_var)
  call write_data(fileobj, "axis_with_edges", data_in_var)
  call write_data(fileobj, "bounds", data_in_var_edges)
  call write_data(fileobj, "edges", data_in_var_edges)

  call close_file(fileobj)

  call open_netcdf_r(fileobj)

  !< Case 1: Here the variable "axis" in the file does not have the attribute "bounds" or "edges", so
  !! it calculates them from the data in "axis"
  answers = 0._k
  call axis_edges(fileobj, "axis", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 1)")

  !< Case 2: Here the variable "axis_with_bounds" in the file has the attribute
  !! "bounds", so the data is read from the variable "bounds"
  answers = 0._k
  call axis_edges(fileobj, "axis_with_bounds", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 2)")

  !< Case 3: Here the variable "axis_with_edges" in the file has the attribute
  !"edges", so the data is read from the variable "edges"
  answers = 0._k
  call axis_edges(fileobj, "axis_with_edges", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 3)")

  !< Case 4: Here the flag "reproduce_null_char_bug_flag" is turned on, so the
  !! edges are calculated from the data in axis because edges has a null character
  !! in the end
  answers = 0._k
  call axis_edges(fileobj, "axis_with_edges", answers, reproduce_null_char_bug_flag=.true.)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 4)")

  call close_file(fileobj)
end subroutine

subroutine test_tranlon
  real(k), dimension(5) :: lon1, lon2, lon3

  lon1 = [1._k, 2._k, 3._k, 4._k,   5._k]
  lon2 = [2._k, 3._k, 4._k, 5._k,   361._k]
  lon3 = [3._k, 4._k, 5._k, 361._k, 362._k]

  ! The first two cases fail due to tranlon's unexpected behavior when no elements are translated.
  ! TODO: Uncomment these tests if/when tranlon's behavior is fixed.

  !call tranlon_assert(lon1, lon1, 0.0_k,    1)
  !call tranlon_assert(lon1, lon1, 1.0_k,    1)

  call tranlon_assert(lon1, lon2, 1.5_k,    2)
  call tranlon_assert(lon1, lon2, 2.0_k,    2)
  call tranlon_assert(lon1, lon3, 2.001_k,  3)
end subroutine

subroutine tranlon_assert(lon0, lon_expected, lon_start, istrt_expected)
  real(k), intent(in) :: lon0(:), lon_expected(:), lon_start
  integer, intent(in) :: istrt_expected
  integer :: istrt_test, i
  real(k) :: lon_test(size(lon0))
  character(:), allocatable :: test_name

  test_name = "tranlon(" // stringify(lon0) // ", " // string(lon_start) // ", istrt)"

  lon_test = lon0
  call tranlon(lon_test, lon_start, istrt_test)
  call array_compare_1d(lon_test, lon_expected, test_name // " unit test failed")

  if (istrt_test.ne.istrt_expected) then
    write(stderr(), "(A)") test_name // " returned erroneous istrt value: " // string(istrt_test)
    write(stderr(), "(A)") "Expected istrt value: " // string(istrt_expected)
    call mpp_error(FATAL, "tranlon unit test failed")
  endif
end subroutine

! Status: SKELETAL
! TODO: More comprehensive interp_1d_1d test
subroutine test_interp_1d_1d
  real(k) :: grid1(8), grid2(5), data1(8), data2(5)

  grid1 = [1._k, 2._k, 3._k, 4._k, 5._k, 6._k, 7._k, 8._k]
  grid2 = [2._k, 3._k, 4._k, 5._k, 6._k]
  data1 = [101._k, 102._k, 103._k, 104._k, 105._k, 106._k, 107._k, 108._k]
  data2 = [102._k, 103._k, 104._k, 105._k, 106._k]

  call interp_1d_1d_assert(grid1, grid2, data1, data2, "linear")
  call interp_1d_1d_assert(grid1, grid2, data1, data2, "cubic_spline")
end subroutine

subroutine interp_1d_1d_assert(grid1, grid2, data1, data2_expected, method, yp1, yp2)
  real(k), intent(in), dimension(:) :: grid1, grid2, data1, data2_expected
  character(*), intent(in), optional :: method
  real(k), intent(in), optional :: yp1, yp2
  real(k) :: data2_test(size(data2_expected))
  character(:), allocatable :: test_name

  test_name = "interp_1d_1d(" // &
              stringify(grid1) // ", " // &
              stringify(grid2) // ", " // &
              stringify(data1) // ", data2"

  if (present(method)) then
    test_name = test_name // ", method=" // method
  endif

  if (present(yp1)) then
    test_name = test_name // ", yp1=" // string(yp1)
  endif

  if (present(yp2)) then
    test_name = test_name // ", yp2=" // string(yp2)
  endif

  test_name = test_name // ")"

  call interp_1d(grid1, grid2, data1, data2_test, method, yp1, yp2)
  call array_compare_1d(data2_test, data2_expected, test_name // " unit test failed")
end subroutine

! Status: SKELETAL
! TODO: More comprehensive interp_1d_2d test
subroutine test_interp_1d_2d
  real(k) :: grid1(2,4), grid2(2,2), data1(2,4), data2(2,2)

  grid1(1,:) = [1._k, 2._k, 3._k, 4._k]
  grid1(2,:) = [5._k, 6._k, 7._k, 8._k]

  grid2(1,:) = [2._k, 3._k]
  grid2(2,:) = [6._k, 7._k]

  data1(1,:) = [101._k, 102._k, 103._k, 104._k]
  data1(2,:) = [105._k, 106._k, 107._k, 108._k]

  data2(1,:) = [102._k, 103._k]
  data2(2,:) = [106._k, 107._k]

  call interp_1d_2d_assert(grid1, grid2, data1, data2)
end subroutine

subroutine interp_1d_2d_assert(grid1, grid2, data1, data2_expected)
  real(k), intent(in), dimension(:,:) :: grid1, grid2, data1, data2_expected
  real(k) :: data2_test(size(data2_expected,1), size(data2_expected,2))
  character(:), allocatable :: test_name

  test_name = "interp_1d_2d(" // &
              stringify(grid1) // ", " // &
              stringify(grid2) // ", " // &
              stringify(data1) // ", data2)"

  call interp_1d(grid1, grid2, data1, data2_test)
  call array_compare_2d(data2_test, data2_expected, test_name // " unit test failed")
end subroutine

! Status: SKELETAL
! TODO: More comprehensive interp_1d_3d test
subroutine test_interp_1d_3d
  real(k) :: grid1(2,2,4), grid2(2,2,2), data1(2,2,4), data2(2,2,2)

  grid1(1,1,:) = [1._k, 2._k, 3._k, 4._k]
  grid1(1,2,:) = [5._k, 6._k, 7._k, 8._k]
  grid1(2,1,:) = [21._k, 22._k, 23._k, 24._k]
  grid1(2,2,:) = [25._k, 26._k, 27._k, 28._k]

  grid2(1,1,:) = [2._k, 3._k]
  grid2(1,2,:) = [6._k, 7._k]
  grid2(2,1,:) = [22._k, 23._k]
  grid2(2,2,:) = [26._k, 27._k]

  data1(1,1,:) = [101._k, 102._k, 103._k, 104._k]
  data1(1,2,:) = [105._k, 106._k, 107._k, 108._k]
  data1(2,1,:) = [201._k, 202._k, 203._k, 204._k]
  data1(2,2,:) = [205._k, 206._k, 207._k, 208._k]

  data2(1,1,:) = [102._k, 103._k]
  data2(1,2,:) = [106._k, 107._k]
  data2(2,1,:) = [202._k, 203._k]
  data2(2,2,:) = [206._k, 207._k]

  call interp_1d_3d_assert(grid1, grid2, data1, data2)
  call interp_1d_3d_assert(grid1, grid2, data1, data2, "linear")
  call interp_1d_3d_assert(grid1, grid2, data1, data2, "cubic_spline")
end subroutine

subroutine interp_1d_3d_assert(grid1, grid2, data1, data2_expected, method, yp1, yp2)
  real(k), intent(in), dimension(:,:,:) :: grid1, grid2, data1, data2_expected
  character(*), intent(in), optional :: method
  real(k), intent(in), optional :: yp1, yp2
  real(k) :: data2_test(size(data2_expected,1), size(data2_expected,2), size(data2_expected,3))
  integer :: i,i2,i3
  character(:), allocatable :: test_name

  test_name = "interp_1d_3d(" // &
              stringify(grid1) // ", " // &
              stringify(grid2) // ", " // &
              stringify(data1) // ", data2"

  if (present(method)) then
    test_name = test_name // ", method=" // method
  endif

  if (present(yp1)) then
    test_name = test_name // ", yp1=" // string(yp1)
  endif

  if (present(yp2)) then
    test_name = test_name // ", yp2=" // string(yp2)
  endif

  test_name = test_name // ")"

  call interp_1d(grid1, grid2, data1, data2_test, method, yp1, yp2)
  call array_compare_3d(data2_test, data2_expected, test_name // " unit test failed")
end subroutine

!
! Supporting utilities
!

subroutine open_netcdf_w(fileobj)
  type(FmsNetcdfFile_t), intent(out) :: fileobj

  if (.not.open_file(fileobj, "test_axis_utils.nc", "overwrite")) then
    call mpp_error(FATAL, "Error opening test_axis_utils.nc to write")
  endif
end subroutine

subroutine open_netcdf_r(fileobj)
  type(FmsNetcdfFile_t), intent(out) :: fileobj

  if (.not.open_file(fileobj, "test_axis_utils.nc", "read")) then
    call mpp_error(FATAL, "Error opening test_axis_utils.nc to read")
  endif
end subroutine

subroutine array_compare_1d(arr1, arr2, msg)
  real(k), intent(in), dimension(:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i, m, n

  m = size(arr1)
  n = size(arr2)

  if (m.ne.n) then
    write(stderr(), "(A)") "1D array comparison failed due to incompatible array sizes"
    write(stderr(), "(A)") "Array 1 has size " // string(m) // " and array 2 has size " // string(n)
    call mpp_error(FATAL, msg)
  endif

  do i=1,m
    if (arr1(i).ne.arr2(i)) then
      write(stderr(), "(A)") "1D array comparison failed due to element " // string(i)
      write(stderr(), "(A)") "Array 1 has value " // string(arr1(i)) // &
                           & " and array 2 has value " // string(arr2(i))
      call mpp_error(FATAL, msg)
    endif
  enddo
end subroutine

subroutine array_compare_2d(arr1, arr2, msg)
  real(k), intent(in), dimension(:,:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i1, i2, m1, m2, n1, n2

  m1 = size(arr1, 1)
  m2 = size(arr1, 2)

  n1 = size(arr2, 1)
  n2 = size(arr2, 2)

  if (m1.ne.n1 .or. m2.ne.n2) then
    write(stderr(), "(A)") "2D array comparison failed due to incompatible array sizes"
    write(stderr(), "(A)") "Array 1 has size " // string(m1) // "x" // string(m2) // &
                          & " and array 2 has size " // string(n1) // "x" // string(n2)
    call mpp_error(FATAL, msg)
  endif

  do i2=1,m2
    do i1=1,m1
      if (arr1(i1,i2).ne.arr2(i1,i2)) then
        write(stderr(), "(A)") "2D array comparison failed due to element " // string(i1) // "," // string(i2)
        write(stderr(), "(A)") "Array 1 has value " // string(arr1(i1,i2)) // &
                             & " and array 2 has value " // string(arr2(i1,i2))
        call mpp_error(FATAL, msg)
      endif
    enddo
  enddo
end subroutine

subroutine array_compare_3d(arr1, arr2, msg)
  real(k), intent(in), dimension(:,:,:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i1, i2, i3, m1, m2, m3, n1, n2, n3

  m1 = size(arr1, 1)
  m2 = size(arr1, 2)
  m3 = size(arr1, 3)

  n1 = size(arr2, 1)
  n2 = size(arr2, 2)
  n3 = size(arr2, 3)

  if (m1.ne.n1 .or. m2.ne.n2 .or. m3.ne.n3) then
    write(stderr(), "(A)") "3D array comparison failed due to incompatible array sizes"
    write(stderr(), "(A)") "Array 1 has size " // string(m1) // "x" // string(m2) // "x" // string(m3) // &
                           & " and array 2 has size " // string(n1) // "x" // string(n2) // "x" // string(n3)
    call mpp_error(FATAL, msg)
  endif

  do i3=1,m3
    do i2=1,m2
      do i1=1,m1
        if (arr1(i1,i2,i3).ne.arr2(i1,i2,i3)) then
          write(stderr(), "(A)") "3D array comparison failed due to element " // &
                               & string(i1) // "," // string(i2) // "," // string(i3)
          write(stderr(), "(A)") "Array 1 has value " // string(arr1(i1,i2,i3)) // &
                               & " and array 2 has value " // string(arr2(i1,i2,i3))
          call mpp_error(FATAL, msg)
        endif
      enddo
    enddo
  enddo
end subroutine

end program test_axis_utils
