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

! Status values:
! * (3/14) TODO     : Not yet implemented
! * (8/14) SKELETAL : Skeletal test has been implemented; comprehensive test has not yet been implemented
! * (3/14) DONE     : Comprehensive test has been implemented

#define PRETTY(x) trim(adjustl(string(x)))

program test_axis_utils

use fms_mod,         only : fms_init, fms_end
use fms2_io_mod, only: FmsNetcdfFile_t, open_file, close_file, register_axis, register_field, &
                     & register_variable_attribute, write_data
use platform_mod, only: r4_kind, r8_kind
use mpp_mod, only: mpp_error, fatal, mpp_pe, mpp_root_pe, mpp_npes, mpp_get_current_pelist, mpp_sync, stdout, stderr
use axis_utils2_mod
use fms_string_utils_mod, only: string

implicit none

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
      write(stderr(),"(2A)") "Unrecognized command line option: ", trim(arg)
  end select
enddo

call fms_end

contains

!
! The actual unit tests
!

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

! Status: TODO
! subroutine get_axis_cart(fileobj, axisname, cart)
subroutine test_get_axis_cart
  type(FmsNetcdfFile_t) :: fileobj

  write(stderr(), "(A)") "Warning: get_axis_cart unit test not yet implemented"
end subroutine

! Status: DONE
subroutine test_lon_in_range
  ! Test some cases where no translation is needed
  call lon_in_range_assert(real(0.0, KIND),          real(0.0, KIND),  real(0.0, KIND))
  call lon_in_range_assert(real(1.0, KIND),          real(0.0, KIND),  real(1.0, KIND))
  call lon_in_range_assert(real(350.0, KIND),        real(0.0, KIND),  real(350.0, KIND))
  call lon_in_range_assert(real(1.0, KIND),          real(1.0, KIND),  real(1.0, KIND))
  call lon_in_range_assert(real(350.0, KIND),        real(1.0, KIND),  real(350.0, KIND))
  call lon_in_range_assert(real(359.0, KIND),        real(0.0, KIND),  real(359.0, KIND))
  call lon_in_range_assert(real(359.0, KIND),        real(1.0, KIND),  real(359.0, KIND))

  ! Test up-translation
  call lon_in_range_assert(real(-2e0, KIND),         real(-1e0, KIND), real(358e0, KIND))
  call lon_in_range_assert(real(-2e0, KIND),         real(0e0, KIND),  real(358e0, KIND))
  call lon_in_range_assert(real(-2e0, KIND),         real(5e0, KIND),  real(358e0, KIND))
  call lon_in_range_assert(real(-1e0, KIND),         real(0e0, KIND),  real(359e0, KIND))
  call lon_in_range_assert(real(-1e0, KIND),         real(5e0, KIND),  real(359e0, KIND))
  call lon_in_range_assert(real(0e0, KIND),          real(5e0, KIND),  real(360e0, KIND))
  call lon_in_range_assert(real(1e0, KIND),          real(5e0, KIND),  real(361e0, KIND))

  ! Test down-translation
  call lon_in_range_assert(real(359.0, KIND),        real(-1.0, KIND), real(-1.0, KIND))
  call lon_in_range_assert(real(360.0, KIND),        real(-1.0, KIND), real(0.0, KIND))
  call lon_in_range_assert(real(360.0, KIND),        real(0.0, KIND),  real(0.0, KIND))
  call lon_in_range_assert(real(361.0, KIND),        real(-1.0, KIND), real(1.0, KIND))
  call lon_in_range_assert(real(361.0, KIND),        real(0.0, KIND),  real(1.0, KIND))
  call lon_in_range_assert(real(362.0, KIND),        real(-1.0, KIND), real(2.0, KIND))
  call lon_in_range_assert(real(362.0, KIND),        real(0.0, KIND),  real(2.0, KIND))

  ! Test rounding behavior
  call lon_in_range_assert(real(1e-5, KIND),         real(0.0, KIND),  real(0.0, KIND))
  call lon_in_range_assert(real(1e-3, KIND),         real(0.0, KIND),  real(1e-3, KIND))
  call lon_in_range_assert(real(360.0 - 1e-5, KIND), real(0.0, KIND),  real(0.0, KIND))
  call lon_in_range_assert(real(360.0 - 1e-3, KIND), real(0.0, KIND),  real(360 - 1e-3, KIND))
end subroutine

subroutine lon_in_range_assert(lon, l_start, ret_expected)
  real(KIND), intent(in) :: lon, l_start, ret_expected
  real(KIND) :: ret_test

  ret_test = lon_in_range(lon, l_start)

  if (ret_test /= ret_expected) then
    write(stderr(), "(6A)") "lon_in_range(", PRETTY(lon), ", ", PRETTY(l_start), ") returned erroneous value: ", PRETTY(ret_test)
    write(stderr(), "(2A)") "Expected return value: ", PRETTY(ret_expected)
    call mpp_error(FATAL, "lon_in_range unit test failed")
  endif
end subroutine

! Status: DONE
subroutine test_frac_index
  real(KIND) :: values(6), v
  integer :: i, n

  values = [1.0, 2.0, 3.0, 5.0, 10.0, 11.0]
  n = size(values)

  ! Test values outside of the input array
  call frac_index_assert(real(values(1) - 0.5, KIND), values, real(-1.0, KIND))
  call frac_index_assert(real(values(n) + 0.5, KIND), values, real(-1.0, KIND))

  ! Test the actual indices
  do i=1,n
    v = values(i)
    call frac_index_assert(v, values, real(i, KIND))
  enddo

  ! Test the 10% point
  do i=1,n-1
    v = values(i) + 0.1*(values(i+1) - values(i))
    call frac_index_assert(v, values, real(i + 0.1, KIND))
  enddo

  ! Test the 25% point
  do i=1,n-1
    v = values(i) + 0.25*(values(i+1) - values(i))
    call frac_index_assert(v, values, real(i + 0.25, KIND))
  enddo

  ! Test the mid-point
  do i=1,n-1
    v = values(i) + 0.5*(values(i+1) - values(i))
    call frac_index_assert(v, values, real(i + 0.5, KIND))
  enddo

  ! Test the 99% point
  do i=1,n-1
    v = values(i) + 0.99*(values(i+1) - values(i))
    call frac_index_assert(v, values, real(i + 0.99, KIND))
  enddo
end subroutine

subroutine frac_index_assert(fval, arr, ret_expected)
  real(KIND), intent(in) :: fval, arr(:), ret_expected
  real(KIND) :: ret_test

  ret_test = frac_index(fval, arr)

  if (ret_test /= ret_expected) then
    write(stderr(), "(4A)") "frac_index(", PRETTY(fval), ", ...) returned erroneous value: ", PRETTY(ret_test)
    write(stderr(), "(2A)") "Expected return value: ", PRETTY(ret_expected)
    call mpp_error(FATAL, "frac_index unit test failed")
  endif
end subroutine

! Status: SKELETAL
subroutine test_frac_index_fail
  real(KIND) :: values(5)
  real(KIND) :: ret_test

  values = [1.0, 2.0, 4.0, 3.0, 5.0]
  ret_test = frac_index(real(1.5, KIND), values)
end subroutine

! Status: SKELETAL
subroutine test_nearest_index
  real(KIND) :: arr(5)

  arr = [5.0, 12.0, 20.0, 40.0, 100.0]

  ! Test values beyond array boundaries
  call nearest_index_assert(real(4.0, KIND), arr, 1)
  call nearest_index_assert(real(1000.0, KIND), arr, size(arr))

  ! Test values actually in the array
  call nearest_index_assert(real(5.0, KIND), arr, 1)
  call nearest_index_assert(real(12.0, KIND), arr, 2)
  call nearest_index_assert(real(20.0, KIND), arr, 3)
  call nearest_index_assert(real(40.0, KIND), arr, 4)
  call nearest_index_assert(real(100.0, KIND), arr, 5)

  ! Test the intervals between array values
  call nearest_index_assert(real(6.0, KIND), arr, 1)
  call nearest_index_assert(real(11.0, KIND), arr, 2)
  call nearest_index_assert(real(15.0, KIND), arr, 2)
  call nearest_index_assert(real(18.0, KIND), arr, 3)
  call nearest_index_assert(real(29.0, KIND), arr, 3)
end subroutine

subroutine nearest_index_assert(val, arr, ret_expected)
  real(KIND), intent(in) :: val, arr(:)
  integer, intent(in) :: ret_expected
  integer :: ret_test

  ret_test = nearest_index(val, arr)

  if (ret_test /= ret_expected) then
    write(stderr(), "(4A)") "nearest_index(", PRETTY(val), ", ...) returned erroneous value: ", PRETTY(ret_test)
    write(stderr(), "(2A)") "Expected return value: ", PRETTY(ret_expected)
    call mpp_error(FATAL, "nearest_index unit test failed")
  endif
end subroutine

! Status: SKELETAL
subroutine test_nearest_index_fail
  real(KIND) :: arr(5)
  integer :: ret_test

  arr=[5.0, 12.0, 40.0, 20.0, 100.0]
  ret_test = nearest_index(real(5.0, KIND), arr)
end subroutine

! Status: DONE
subroutine test_axis_edges
  real(KIND) :: data_in_var(10)
  real(KIND) :: data_in_var_edges(2,10)
  real(KIND) :: data_in_answers(11)
  type(FmsNetcdfFile_t) :: fileobj
  real(KIND)    :: answers(11)
  integer :: i

  do i=1,10
     data_in_var(i) = real(i, KIND) - 0.5

     data_in_var_edges(1,i) = real(i-1, KIND)
     data_in_var_edges(2,i) = real(i, KIND)

     data_in_answers(i) = real(i-1, KIND)
  enddo

  data_in_answers(11) = real(10, KIND)

  if (mpp_pe() .eq. mpp_root_pe()) then
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
  endif

  call mpp_sync

  call open_netcdf_r(fileobj)

  !< Case 1: Here the variable "axis" in the file does not have the attribute "bounds" or "edges", so
  !! it calculates them from the data in "axis"
  answers = 0.0
  call axis_edges(fileobj, "axis", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 1)")

  !< Case 2: Here the variable "axis_with_bounds" in the file has the attribute
  !! "bounds", so the data is read from the variable "bounds"
  answers = 0.0
  call axis_edges(fileobj, "axis_with_bounds", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 2)")

  !< Case 3: Here the variable "axis_with_edges" in the file has the attribute
  !"edges", so the data is read from the variable "edges"
  answers = 0.0
  call axis_edges(fileobj, "axis_with_edges", answers)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 3)")

  !< Case 4: Here the flag "reproduce_null_char_bug_flag" is turned on, so the
  !! edges are calculated from the data in axis because edges has a null character
  !! in the end
  answers = 0.0
  call axis_edges(fileobj, "axis_with_edges", answers, reproduce_null_char_bug_flag=.true.)
  call array_compare_1d(answers, data_in_answers, "axis_edges unit test failed (case 4)")

  call close_file(fileobj)
end subroutine

! Status: SKELETAL
subroutine test_tranlon
  real(KIND), dimension(5) :: lon1, lon2, lon3

  lon1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  lon2 = [2.0, 3.0, 4.0, 5.0, 361.0]
  lon3 = [3.0, 4.0, 5.0, 361.0, 362.0]

  ! The first two cases seem to reveal an error in tranlon. Should tranlon be changed so that
  ! istrt=1 in the first two cases?
  call tranlon_assert(lon1, lon1, real(0.0, KIND), 0)
  call tranlon_assert(lon1, lon1, real(1.0, KIND), 0)
  call tranlon_assert(lon1, lon2, real(1.5, KIND), 2)
  call tranlon_assert(lon1, lon2, real(2.0, KIND), 2)
  call tranlon_assert(lon1, lon3, real(2.001, KIND), 3)
end subroutine

subroutine tranlon_assert(lon0, lon_expected, lon_start, istrt_expected)
  real(KIND), intent(in) :: lon0(:), lon_expected(:), lon_start
  integer, intent(in) :: istrt_expected
  integer :: istrt_test, i
  real(KIND) :: lon_test(size(lon0))

  lon_test = lon0
  call tranlon(lon_test, lon_start, istrt_test)
  call array_compare_1d(lon_test, lon_expected, "tranlon unit test failed")

  if (istrt_test.ne.istrt_expected) then
    write(stderr(), "(4A)") "tranlon(..., ", PRETTY(lon_start), ", istrt) returned erroneous istrt value: ", PRETTY(istrt_test)
    write(stderr(), "(2A)") "Expected istrt value: ", PRETTY(istrt_expected)
    call mpp_error(FATAL, "tranlon unit test failed")
  endif
end subroutine

! Status: SKELETAL
! subroutine interp_1d_1d(grid1,grid2,data1,data2, method, yp1, yp2)
subroutine test_interp_1d_1d
  real(KIND) :: grid1(8), grid2(5), data1(8), data2(5)

  grid1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  grid2 = [2.0, 3.0, 4.0, 5.0, 6.0]
  data1 = [101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0]
  data2 = [102.0, 103.0, 104.0, 105.0, 106.0]

  call interp_1d_1d_assert(grid1, grid2, data1, data2, "linear")
  call interp_1d_1d_assert(grid1, grid2, data1, data2, "cubic_spline")
end subroutine

subroutine interp_1d_1d_assert(grid1, grid2, data1, data2_expected, method, yp1, yp2)
  real(KIND), intent(in), dimension(:) :: grid1, grid2, data1, data2_expected
  character(*), intent(in), optional :: method
  real(KIND), intent(in), optional :: yp1, yp2
  real(KIND) :: data2_test(size(data2_expected))

  call interp_1d(grid1, grid2, data1, data2_test, method, yp1, yp2)
  call array_compare_1d(data2_test, data2_expected, "interp_1d_1d unit test failed")
end subroutine

! Status: SKELETAL
subroutine test_interp_1d_2d
  real(KIND) :: grid1(2,4), grid2(2,2), data1(2,4), data2(2,2)

  grid1(1,:) = [1.0, 2.0, 3.0, 4.0]
  grid1(2,:) = [5.0, 6.0, 7.0, 8.0]

  grid2(1,:) = [2.0, 3.0]
  grid2(2,:) = [6.0, 7.0]

  data1(1,:) = [101.0, 102.0, 103.0, 104.0]
  data1(2,:) = [105.0, 106.0, 107.0, 108.0]

  data2(1,:) = [102.0, 103.0]
  data2(2,:) = [106.0, 107.0]

  call interp_1d_2d_assert(grid1, grid2, data1, data2)
end subroutine

subroutine interp_1d_2d_assert(grid1, grid2, data1, data2_expected)
  real(KIND), intent(in), dimension(:,:) :: grid1, grid2, data1, data2_expected
  real(KIND) :: data2_test(size(data2_expected,1), size(data2_expected,2))

  call interp_1d(grid1, grid2, data1, data2_test)
  call array_compare_2d(data2_test, data2_expected, "interp_1d_2d unit test failed")
end subroutine

! Status: SKELETAL
subroutine test_interp_1d_3d
  real(KIND) :: grid1(2,2,4), grid2(2,2,2), data1(2,2,4), data2(2,2,2)

  grid1(1,1,:) = [1.0, 2.0, 3.0, 4.0]
  grid1(1,2,:) = [5.0, 6.0, 7.0, 8.0]
  grid1(2,1,:) = [21.0, 22.0, 23.0, 24.0]
  grid1(2,2,:) = [25.0, 26.0, 27.0, 28.0]

  grid2(1,1,:) = [2.0, 3.0]
  grid2(1,2,:) = [6.0, 7.0]
  grid2(2,1,:) = [22.0, 23.0]
  grid2(2,2,:) = [26.0, 27.0]

  data1(1,1,:) = [101.0, 102.0, 103.0, 104.0]
  data1(1,2,:) = [105.0, 106.0, 107.0, 108.0]
  data1(2,1,:) = [201.0, 202.0, 203.0, 204.0]
  data1(2,2,:) = [205.0, 206.0, 207.0, 208.0]

  data2(1,1,:) = [102.0, 103.0]
  data2(1,2,:) = [106.0, 107.0]
  data2(2,1,:) = [202.0, 203.0]
  data2(2,2,:) = [206.0, 207.0]

  call interp_1d_3d_assert(grid1, grid2, data1, data2)
  call interp_1d_3d_assert(grid1, grid2, data1, data2, "linear")
  call interp_1d_3d_assert(grid1, grid2, data1, data2, "cubic_spline")
end subroutine

subroutine interp_1d_3d_assert(grid1, grid2, data1, data2_expected, method, yp1, yp2)
  real(KIND), intent(in), dimension(:,:,:) :: grid1, grid2, data1, data2_expected
  character(*), intent(in), optional :: method
  real(KIND), intent(in), optional :: yp1, yp2
  real(KIND) :: data2_test(size(data2_expected,1), size(data2_expected,2), size(data2_expected,3))
  integer :: i,j,k

  call interp_1d(grid1, grid2, data1, data2_test, method, yp1, yp2)
  call array_compare_3d(data2_test, data2_expected, "interp_1d_3d unit test failed")
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
  integer, allocatable  :: pes(:)

  allocate(pes(mpp_npes()))
  call mpp_get_current_pelist(pes)

  if (.not.open_file(fileobj, "test_axis_utils.nc", "read", pelist=pes)) then
    call mpp_error(FATAL, "Error opening test_axis_utils.nc to read")
  endif

  deallocate(pes)
end subroutine

subroutine array_compare_1d(arr1, arr2, msg)
  real(KIND), intent(in), dimension(:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i, n, n2

  n = size(arr1)
  n2 = size(arr2)

  if (n2.ne.n) then
    write(stderr(), "(A)") "1D array comparison failed due to incompatible array sizes"
    write(stderr(), "(4A)") "Array 1 has size ", PRETTY(n), " and array 2 has size ", PRETTY(n2)
    call mpp_error(FATAL, msg)
  endif

  do i=1,n
    if (arr1(i).ne.arr2(i)) then
      write(stderr(), "(2A)") "1D array comparison failed due to element ", PRETTY(i)
      write(stderr(), "(4A)") "Array 1 has value ", PRETTY(arr1(i)), " and array 2 has value ", PRETTY(arr2(i))
      call mpp_error(FATAL, msg)
    endif
  enddo
end subroutine

subroutine array_compare_2d(arr1, arr2, msg)
  real(KIND), intent(in), dimension(:,:) :: arr1, arr2
  character(*), intent(in) :: msg
  integer :: i,j,m,n,m2,n2

  m = size(arr1, 1)
  n = size(arr1, 2)

  m2 = size(arr2, 1)
  n2 = size(arr2, 2)

  if (m.ne.m2 .or. n.ne.n2) then
    write(stderr(), "(A)") "2D array comparison failed due to incompatible array sizes"
    write(stderr(), "(8A)") "Array 1 has size ", PRETTY(m), "x", PRETTY(n), &
                          & " and array 2 has size ", PRETTY(m2), "x", PRETTY(n2)
    call mpp_error(FATAL, msg)
  endif

  do i=1,n
    do j=1,m
      if (arr1(j,i).ne.arr2(j,i)) then
        write(stderr(), "(4A)") "2D array comparison failed due to element ", PRETTY(j), ",", PRETTY(i)
        write(stderr(), "(4A)") "Array 1 has value ", PRETTY(arr1(j,i)), " and array 2 has value ", PRETTY(arr2(j,i))
        call mpp_error(FATAL, msg)
      endif
    enddo
  enddo
end subroutine

subroutine array_compare_3d(arr1, arr2, msg)
  real(KIND), intent(in), dimension(:,:,:) :: arr1, arr2
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
    write(stderr(), "(12A)") "Array 1 has size ", PRETTY(l), "x", PRETTY(m), "x", PRETTY(n), &
                           & " and array 2 has size ", PRETTY(l2), "x", PRETTY(m2), "x", PRETTY(n2)
    call mpp_error(FATAL, msg)
  endif

  do i=1,n
    do j=1,m
      do k=1,l
        if (arr1(k,j,i).ne.arr2(k,j,i)) then
          write(stderr(), "(6A)") "3D array comparison failed due to element ", PRETTY(k), ",", PRETTY(j), ",", PRETTY(i)
          write(stderr(), "(4A)") "Array 1 has value ", PRETTY(arr1(k,j,i)), " and array 2 has value ", PRETTY(arr2(k,j,i))
          call mpp_error(FATAL, msg)
        endif
      enddo
    enddo
  enddo
end subroutine

end program test_axis_utils
