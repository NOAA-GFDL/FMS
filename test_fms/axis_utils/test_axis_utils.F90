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

! gfortran lacks support for the macro pasting operator, but it does support
! whitespace around the underscore.
#ifdef __GFORTRAN__
#define C(x) x _ AU_TEST_KIND_
#else
#define C(x) x ## _ ## AU_TEST_KIND_
#endif

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

    case ('--nearest-index')
      print "(A)", "Testing nearest_index"
      call test_nearest_index

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

    ! Check an unknown value on a "cartesian_axis" or "axis" attribute
    ! TODO: This test fails. Should get_axis_cart be changed, or should this test be changed?
    attr_value = "unexpected"
    var_name = attr_name // "_attr_value_" // attr_value
    call get_axis_cart_test_add(test, var_name, "N")
    call register_variable_attribute(test%fileobj, var_name, attr_name, attr_value, str_len=len(attr_value))

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

#define r4_kind "float"
#define r8_kind "double"
  character(*), parameter :: kind_str = AU_TEST_KIND_
#undef r4_kind
#undef r8_kind

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
  real(AU_TEST_KIND_), parameter :: eps_big = C(1e-3), eps_tiny = C(1e-5)

  ! Test some cases where no translation is needed
  call lon_in_range_assert(C(0.),              C(0.),  C(0.))
  call lon_in_range_assert(C(1.),              C(0.),  C(1.))
  call lon_in_range_assert(C(350.),            C(0.),  C(350.))
  call lon_in_range_assert(C(1.),              C(1.),  C(1.))
  call lon_in_range_assert(C(350.),            C(1.),  C(350.))
  call lon_in_range_assert(C(359.),            C(0.),  C(359.))
  call lon_in_range_assert(C(359.),            C(1.),  C(359.))

  ! Test up-translation
  call lon_in_range_assert(C(-2.),             C(-1.), C(358.))
  call lon_in_range_assert(C(-2.),             C(0.),  C(358.))
  call lon_in_range_assert(C(-2.),             C(5.),  C(358.))
  call lon_in_range_assert(C(-1.),             C(0.),  C(359.))
  call lon_in_range_assert(C(-1.),             C(5.),  C(359.))
  call lon_in_range_assert(C(0.),              C(5.),  C(360.))
  call lon_in_range_assert(C(1.),              C(5.),  C(361.))

  ! Test down-translation
  call lon_in_range_assert(C(359.),            C(-1.), C(-1.))
  call lon_in_range_assert(C(360.),            C(-1.), C(0.))
  call lon_in_range_assert(C(360.),            C(0.),  C(0.))
  call lon_in_range_assert(C(361.),            C(-1.), C(1.))
  call lon_in_range_assert(C(361.),            C(0.),  C(1.))
  call lon_in_range_assert(C(362.),            C(-1.), C(2.))
  call lon_in_range_assert(C(362.),            C(0.),  C(2.))

  ! Test rounding behavior
  call lon_in_range_assert(eps_tiny,           C(0.),  C(0.))
  call lon_in_range_assert(eps_big,            C(0.),  eps_big)
  call lon_in_range_assert(C(360.) - eps_tiny, C(0.),  C(0.))
  call lon_in_range_assert(C(360.) - eps_big,  C(0.),  C(360.) - eps_big)
end subroutine

subroutine lon_in_range_assert(lon, l_start, ret_expected)
  real(AU_TEST_KIND_), intent(in) :: lon, l_start, ret_expected
  real(AU_TEST_KIND_) :: ret_test

  ret_test = lon_in_range(lon, l_start)

  if (ret_test /= ret_expected) then
    write(stderr(), "(A)") "lon_in_range(" // string(lon) // ", " // string(l_start) // &
                         & ") returned erroneous value: " // string(ret_test)
    write(stderr(), "(A)") "Expected return value: " // string(ret_expected)
    call mpp_error(FATAL, "lon_in_range unit test failed")
  endif
end subroutine

#define CALC_FRAC_INDEX(i, v, values) real(i, AU_TEST_KIND_) + (v - values(i)) / (values(i + 1) - values(i))

subroutine test_frac_index
  real(AU_TEST_KIND_) :: values(6), v, fi
  integer :: i, n
  real(AU_TEST_KIND_), parameter :: f10=C(0.1), f25=C(0.25), f50=C(0.5), f99=C(0.99)

  values = [C(1.), C(2.), C(3.), C(5.), C(10.), C(11.)]
  n = size(values)

  ! Test values outside of the input array
  call frac_index_assert(real(values(1), AU_TEST_KIND_) - f50, values, C(-1.))
  call frac_index_assert(real(values(n), AU_TEST_KIND_) + f50, values, C(-1.))

  ! Test the actual indices
  do i=1,n
    v = values(i)
    call frac_index_assert(v, values, real(i, AU_TEST_KIND_))
  enddo

  ! Test the 10% point
  do i=1,n-1
    v = values(i) + f10*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo

  ! Test the 25% point
  do i=1,n-1
    v = values(i) + f25*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo

  ! Test the mid-point
  do i=1,n-1
    v = values(i) + f50*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo

  ! Test the 99% point
  do i=1,n-1
    v = values(i) + f99*(values(i+1) - values(i))
    fi = CALC_FRAC_INDEX(i, v, values)
    call frac_index_assert(v, values, fi)
  enddo
end subroutine

subroutine frac_index_assert(fval, arr, ret_expected)
  real(AU_TEST_KIND_), intent(in) :: fval, arr(:), ret_expected
  real(AU_TEST_KIND_) :: ret_test

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
  real(AU_TEST_KIND_) :: values(5)
  real(AU_TEST_KIND_) :: ret_test

  values = [C(1.), C(2.), C(4.), C(3.), C(5.)]
  ret_test = frac_index(C(1.5), values)
end subroutine

subroutine test_nearest_index
  real(AU_TEST_KIND_) :: arr(5)

  arr = [C(5.), C(12.), C(20.), C(40.), C(100.)]

  ! Test values beyond array boundaries
  call nearest_index_assert(C(4.), arr, 1)
  call nearest_index_assert(C(1000.), arr, size(arr))

  ! Test values actually in the array
  call nearest_index_assert(C(5.), arr, 1)
  call nearest_index_assert(C(12.), arr, 2)
  call nearest_index_assert(C(20.), arr, 3)
  call nearest_index_assert(C(40.), arr, 4)
  call nearest_index_assert(C(100.), arr, 5)

  ! Test the intervals between array values
  call nearest_index_assert(C(6.), arr, 1)
  call nearest_index_assert(C(11.), arr, 2)
  call nearest_index_assert(C(15.), arr, 2)
  call nearest_index_assert(C(18.), arr, 3)
  call nearest_index_assert(C(29.), arr, 3)
end subroutine

subroutine nearest_index_assert(val, arr, ret_expected)
  real(AU_TEST_KIND_), intent(in) :: val, arr(:)
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
  real(AU_TEST_KIND_) :: arr(5)
  integer :: ret_test

  arr=[C(5.), C(12.), C(40.), C(20.), C(100.)]
  ret_test = nearest_index(C(5.), arr)
end subroutine

subroutine test_axis_edges
  real(AU_TEST_KIND_) :: data_in_var(10)
  real(AU_TEST_KIND_) :: data_in_var_edges(2,10)
  real(AU_TEST_KIND_) :: data_in_answers(11)
  type(FmsNetcdfFile_t) :: fileobj
  real(AU_TEST_KIND_)    :: answers(11)
  integer :: i

  do i=1,10
     data_in_var(i) = real(i, AU_TEST_KIND_) - C(0.5)

     data_in_var_edges(1,i) = real(i-1, AU_TEST_KIND_)
     data_in_var_edges(2,i) = real(i, AU_TEST_KIND_)

     data_in_answers(i) = real(i-1, AU_TEST_KIND_)
  enddo

  data_in_answers(11) = C(10.)

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
  answers = C(0.)
  call axis_edges(fileobj, "axis", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 1)")

  !< Case 2: Here the variable "axis_with_bounds" in the file has the attribute
  !! "bounds", so the data is read from the variable "bounds"
  answers = C(0.)
  call axis_edges(fileobj, "axis_with_bounds", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 2)")

  !< Case 3: Here the variable "axis_with_edges" in the file has the attribute
  !"edges", so the data is read from the variable "edges"
  answers = C(0.)
  call axis_edges(fileobj, "axis_with_edges", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 3)")

  !< Case 4: Here the flag "reproduce_null_char_bug_flag" is turned on, so the
  !! edges are calculated from the data in axis because edges has a null character
  !! in the end
  answers = C(0.)
  call axis_edges(fileobj, "axis_with_edges", answers, reproduce_null_char_bug_flag=.true.)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 4)")

  call close_file(fileobj)
end subroutine

subroutine test_tranlon
  real(AU_TEST_KIND_), dimension(5) :: lon1, lon2, lon3

  lon1 = [C(1.), C(2.), C(3.), C(4.), C(5.)]
  lon2 = [C(2.), C(3.), C(4.), C(5.), C(361.)]
  lon3 = [C(3.), C(4.), C(5.), C(361.), C(362.)]

  ! TODO: The first two cases fail due to tranlon's unexpected behavior when no elements are translated.
  ! Should tranlon be changed so that istrt=1 in the first two cases, or should the test be changed?
  call tranlon_assert(lon1, lon1, C(0.0),    1)
  call tranlon_assert(lon1, lon1, C(1.0),    1)
  call tranlon_assert(lon1, lon2, C(1.5),    2)
  call tranlon_assert(lon1, lon2, C(2.0),    2)
  call tranlon_assert(lon1, lon3, C(2.001),  3)
end subroutine

subroutine tranlon_assert(lon0, lon_expected, lon_start, istrt_expected)
  real(AU_TEST_KIND_), intent(in) :: lon0(:), lon_expected(:), lon_start
  integer, intent(in) :: istrt_expected
  integer :: istrt_test, i
  real(AU_TEST_KIND_) :: lon_test(size(lon0))
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
  real(AU_TEST_KIND_) :: grid1(8), grid2(5), data1(8), data2(5)

  grid1 = [C(1.), C(2.), C(3.), C(4.), C(5.), C(6.), C(7.), C(8.)]
  grid2 = [C(2.), C(3.), C(4.), C(5.), C(6.)]
  data1 = [C(101.), C(102.), C(103.), C(104.), C(105.), C(106.), C(107.), C(108.)]
  data2 = [C(102.), C(103.), C(104.), C(105.), C(106.)]

  call interp_1d_1d_assert(grid1, grid2, data1, data2, "linear")
  call interp_1d_1d_assert(grid1, grid2, data1, data2, "cubic_spline")
end subroutine

subroutine interp_1d_1d_assert(grid1, grid2, data1, data2_expected, method, yp1, yp2)
  real(AU_TEST_KIND_), intent(in), dimension(:) :: grid1, grid2, data1, data2_expected
  character(*), intent(in), optional :: method
  real(AU_TEST_KIND_), intent(in), optional :: yp1, yp2
  real(AU_TEST_KIND_) :: data2_test(size(data2_expected))
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
  real(AU_TEST_KIND_) :: grid1(2,4), grid2(2,2), data1(2,4), data2(2,2)

  grid1(1,:) = [C(1.), C(2.), C(3.), C(4.)]
  grid1(2,:) = [C(5.), C(6.), C(7.), C(8.)]

  grid2(1,:) = [C(2.), C(3.)]
  grid2(2,:) = [C(6.), C(7.)]

  data1(1,:) = [C(101.), C(102.), C(103.), C(104.)]
  data1(2,:) = [C(105.), C(106.), C(107.), C(108.)]

  data2(1,:) = [C(102.), C(103.)]
  data2(2,:) = [C(106.), C(107.)]

  call interp_1d_2d_assert(grid1, grid2, data1, data2)
end subroutine

subroutine interp_1d_2d_assert(grid1, grid2, data1, data2_expected)
  real(AU_TEST_KIND_), intent(in), dimension(:,:) :: grid1, grid2, data1, data2_expected
  real(AU_TEST_KIND_) :: data2_test(size(data2_expected,1), size(data2_expected,2))
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
  real(AU_TEST_KIND_) :: grid1(2,2,4), grid2(2,2,2), data1(2,2,4), data2(2,2,2)

  grid1(1,1,:) = [C(1.), C(2.), C(3.), C(4.)]
  grid1(1,2,:) = [C(5.), C(6.), C(7.), C(8.)]
  grid1(2,1,:) = [C(21.), C(22.), C(23.), C(24.)]
  grid1(2,2,:) = [C(25.), C(26.), C(27.), C(28.)]

  grid2(1,1,:) = [C(2.), C(3.)]
  grid2(1,2,:) = [C(6.), C(7.)]
  grid2(2,1,:) = [C(22.), C(23.)]
  grid2(2,2,:) = [C(26.), C(27.)]

  data1(1,1,:) = [C(101.), C(102.), C(103.), C(104.)]
  data1(1,2,:) = [C(105.), C(106.), C(107.), C(108.)]
  data1(2,1,:) = [C(201.), C(202.), C(203.), C(204.)]
  data1(2,2,:) = [C(205.), C(206.), C(207.), C(208.)]

  data2(1,1,:) = [C(102.), C(103.)]
  data2(1,2,:) = [C(106.), C(107.)]
  data2(2,1,:) = [C(202.), C(203.)]
  data2(2,2,:) = [C(206.), C(207.)]

  call interp_1d_3d_assert(grid1, grid2, data1, data2)
  call interp_1d_3d_assert(grid1, grid2, data1, data2, "linear")
  call interp_1d_3d_assert(grid1, grid2, data1, data2, "cubic_spline")
end subroutine

subroutine interp_1d_3d_assert(grid1, grid2, data1, data2_expected, method, yp1, yp2)
  real(AU_TEST_KIND_), intent(in), dimension(:,:,:) :: grid1, grid2, data1, data2_expected
  character(*), intent(in), optional :: method
  real(AU_TEST_KIND_), intent(in), optional :: yp1, yp2
  real(AU_TEST_KIND_) :: data2_test(size(data2_expected,1), size(data2_expected,2), size(data2_expected,3))
  integer :: i,j,k
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
  real(AU_TEST_KIND_), intent(in), dimension(:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i, n, n2

  n = size(arr1)
  n2 = size(arr2)

  if (n2.ne.n) then
    write(stderr(), "(A)") "1D array comparison failed due to incompatible array sizes"
    write(stderr(), "(A)") "Array 1 has size " // string(n) // " and array 2 has size " // string(n2)
    call mpp_error(FATAL, msg)
  endif

  do i=1,n
    if (arr1(i).ne.arr2(i)) then
      write(stderr(), "(A)") "1D array comparison failed due to element " // string(i)
      write(stderr(), "(A)") "Array 1 has value " // string(arr1(i)) // &
                           & " and array 2 has value " // string(arr2(i))
      call mpp_error(FATAL, msg)
    endif
  enddo
end subroutine

subroutine array_compare_2d(arr1, arr2, msg)
  real(AU_TEST_KIND_), intent(in), dimension(:,:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i,j,m,n,m2,n2

  m = size(arr1, 1)
  n = size(arr1, 2)

  m2 = size(arr2, 1)
  n2 = size(arr2, 2)

  if (m.ne.m2 .or. n.ne.n2) then
    write(stderr(), "(A)") "2D array comparison failed due to incompatible array sizes"
    write(stderr(), "(A)") "Array 1 has size " // string(m) // "x" // string(n) // &
                          & " and array 2 has size " // string(m2) // "x" // string(n2)
    call mpp_error(FATAL, msg)
  endif

  do i=1,n
    do j=1,m
      if (arr1(j,i).ne.arr2(j,i)) then
        write(stderr(), "(A)") "2D array comparison failed due to element " // string(j) // "," // string(i)
        write(stderr(), "(A)") "Array 1 has value " // string(arr1(j,i)) // &
                             & " and array 2 has value " // string(arr2(j,i))
        call mpp_error(FATAL, msg)
      endif
    enddo
  enddo
end subroutine

subroutine array_compare_3d(arr1, arr2, msg)
  real(AU_TEST_KIND_), intent(in), dimension(:,:,:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i,j,k,l,m,n,l2,m2,n2

  l = size(arr1, 1)
  m = size(arr1, 2)
  n = size(arr1, 3)

  l2 = size(arr2, 1)
  m2 = size(arr2, 2)
  n2 = size(arr2, 3)

  if (l.ne.l2 .or. m.ne.m2 .or. n.ne.n2) then
    write(stderr(), "(A)") "3D array comparison failed due to incompatible array sizes"
    write(stderr(), "(A)") "Array 1 has size " // string(l) // "x" // string(m) // "x" // string(n) // &
                           & " and array 2 has size " // string(l2) // "x" // string(m2) // "x" // string(n2)
    call mpp_error(FATAL, msg)
  endif

  do i=1,n
    do j=1,m
      do k=1,l
        if (arr1(k,j,i).ne.arr2(k,j,i)) then
          write(stderr(), "(A)") "3D array comparison failed due to element " // &
                               & string(k) // "," // string(j) // "," // string(i)
          write(stderr(), "(A)") "Array 1 has value " // string(arr1(k,j,i)) // &
                               & " and array 2 has value " // string(arr2(k,j,i))
          call mpp_error(FATAL, msg)
        endif
      enddo
    enddo
  enddo
end subroutine

end program test_axis_utils
