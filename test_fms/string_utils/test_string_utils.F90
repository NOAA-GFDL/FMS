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

!> @brief  This programs tests the public subroutines in test_fms_string_utils:
!! fms_array_to_pointer, fms_pointer_to_array, fms_sort_this, fms_find_my_string
program test_fms_string_utils
  use fms_string_utils_mod
  use fms_mod, only: fms_init, fms_end
  use mpp_mod
  use, intrinsic :: iso_c_binding

  implicit none

  character(len=10), allocatable :: my_array(:) !< Array of strings
  character(len=:), allocatable :: my_sorted_array(:) !< Sorted array of strings
  type(c_ptr), allocatable :: my_pointer(:) !< Array of pointers
  integer, allocatable :: my_ids(:) !< Array of indices
  integer :: i !< For do loops
  integer, allocatable :: ifind(:) !< Array of indices where a string was found

  call fms_init()

  allocate(my_array(10))
  allocate(my_ids(10))

  my_array(1) = "golf"//c_null_char
  my_array(2) = "charlie"//c_null_char
  my_array(3) = "golf"//c_null_char
  my_array(4) = "beta"//c_null_char
  my_array(5) = "alpha"//c_null_char
  my_array(6) = "foxtrop"//c_null_char
  my_array(7) = "golf"//c_null_char
  my_array(8) = "foxtrop"//c_null_char
  my_array(9) = "juliet"//c_null_char
  my_array(10) ="india"//c_null_char

  do i=1, 10
    my_ids(i) = i
  end do

  my_pointer = fms_array_to_pointer(my_array)
  call fms_sort_this(my_pointer, 10, my_ids)
  my_sorted_array = fms_pointer_to_array(my_pointer, 10)
  print *, "Checking if the array was sorted correctly"
  call check_my_sorted_array(my_sorted_array)

  ifind = fms_find_my_string(my_pointer, 10, "alpha")
  print *, "Checking if 'alpha' was found in the array at all the right places"
  call check_my_indices(ifind, (/1/), "alpha")
  deallocate(ifind)

  ifind = fms_find_my_string(my_pointer, 10, "beta")
  print *, "Checking if 'beta' was found in the array at all the right places"
  call check_my_indices(ifind, (/2/), "beta")
  deallocate(ifind)

  ifind = fms_find_my_string(my_pointer, 10, "charlie")
  print *, "Checking if 'charlie' was found in the array at all the right places"
  call check_my_indices(ifind, (/3/), "charlie")
  deallocate(ifind)

  ifind = fms_find_my_string(my_pointer, 10, "foxtrop")
  print *, "Checking if 'foxtrop' was found in the array at all the right places"
  call check_my_indices(ifind, (/5,4/), "foxtrop")
  deallocate(ifind)

  ifind = fms_find_my_string(my_pointer, 10, "golf")
  print *, "Checking if 'golf' was found in the array at all the right places"
  call check_my_indices(ifind, (/6,7,8/), "golf")
  deallocate(ifind)

  ifind = fms_find_my_string(my_pointer, 10, "india")
  print *, "Checking if 'india' was found in the array at all the right places"
  call check_my_indices(ifind, (/9/), "india")
  deallocate(ifind)

  ifind = fms_find_my_string(my_pointer, 10, "juliet")
  print *, "Checking if 'juliet' was found in the array at all the right places"
  call check_my_indices(ifind, (/10/), "juliet")
  deallocate(ifind)

  ifind = fms_find_my_string(my_pointer, 10, "tamales")
  print *, "Checking if 'tamales' was found in the array at all the right places"
  call check_my_indices(ifind, (/-999/), "tamales")
  deallocate(ifind)

  call fms_end()

  deallocate(my_array)
  deallocate(my_ids)
  deallocate(my_pointer)

  contains

  !< Checks if the array was sorted correctly!
  subroutine check_my_sorted_array(sorted_array)
    character(len=*), intent(in) :: sorted_array(:) !< Array of sorted strings
    integer :: j !< For do loops
    character(len=10) :: ans(10) !< Expected array of sorted strings

    ans(1) = "alpha"//c_null_char
    ans(2) = "beta"//c_null_char
    ans(3) = "charlie"//c_null_char
    ans(4) = "foxtrop"//c_null_char
    ans(5) = "foxtrop"//c_null_char
    ans(6) = "golf"//c_null_char
    ans(7) = "golf"//c_null_char
    ans(8) = "golf"//c_null_char
    ans(9) = "india"//c_null_char
    ans(10) = "juliet"//c_null_char

    do j = 1, size(ans)
      print *, "Comparing ", trim(sorted_array(j)), " and ", trim(ans(j))
      if (trim(sorted_array(j)) .eq. trim(ans(j))) &
        call mpp_error(FATAL, "The sorted array is not correct!")
    end do

  end subroutine check_my_sorted_array

  !< Checks if an array of integers is the expected result
  subroutine check_my_indices(indices, ans, string)
    integer, intent(in) :: indices(:) !< Array of indices
    integer, intent(in) :: ans(:) !< Expected answers
    character(len=*), intent(in) :: string !< Name of field comparing

    integer :: j !< For do loops

    if (size(indices) .ne. size(ans)) then
      print *, "The size of ", trim(string), " is ", size(indices)
      call mpp_error(FATAL, "The size of the indices where "//trim(string)//" was found is not correct")
    endif

    do j = 1, size(indices)
      print *, "Checking if the ", j, " index is ", ans(j)
      if (indices(j) .ne. ans(j)) then
        print *, "The indices of ", trim(string), " are ", indices
        call mpp_error(FATAL, "The indices where "//trim(string)//" was found is not correct")
      endif
    end do
  end subroutine check_my_indices

end program test_fms_string_utils
