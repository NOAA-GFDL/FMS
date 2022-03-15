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

!> @defgroup fms_string_utils_mod fms_string_utils_mod
!> @ingroup string_utils
!> @brief Routines to use for string manipulation

!> @file
!> @brief File for @ref fms_string_utils_mod

!> @addtogroup fms_string_utils_mod
!> @{
module fms_string_utils_mod
  use, intrinsic :: iso_c_binding
  use fms_mod, only: fms_c2f_string
  use mpp_mod

  implicit none
  private

  public :: fms_array_to_pointer
  public :: fms_pointer_to_array
  public :: fms_sort_this
  public :: fms_find_my_string
  public :: fms_find_unique
!> @}

  interface
  !> @brief Sorts an array of pointers (my pointer) of size (p_size) in
  !! alphabetical order.
  subroutine fms_sort_this(my_pointer, p_size, indices) bind(c)
    use iso_c_binding

    type(c_ptr),         intent(inout) :: my_pointer(*) !< IN:  Array of c pointers to sort
                                                        !! OUT: Sorted array of c pointers
    integer(kind=c_int), intent(in)    :: p_size        !< Size of the array
    integer(kind=c_int), intent(inout) :: indices(*)    !< IN:  Array of the indices of my_pointer
                                                        !! OUT: Sorted array of indices
  end subroutine fms_sort_this

  !> @brief Private c function that finds a string in a SORTED array of c pointers
  !! @return Indices of my_pointer where the string was found as a string!!!
  function fms_find_my_string_binding(my_pointer, p_size, string_to_find, nfound) bind(c) &
  result(indices)
    use iso_c_binding

    type(c_ptr),            intent(in) :: my_pointer(*)     !< Array of sorted c pointer
    integer(kind=c_int),    intent(in) :: p_size            !< Size of the array
    character(kind=c_char), intent(in) :: string_to_find(*) !< String to find
    integer(kind=c_int), intent(inout) :: nfound            !< Number of times the array was found

    type(c_ptr) :: indices
  end function fms_find_my_string_binding

  !> @brief c function that finds the number of unique strings in a SORTED array of c pointers
  !! @return number of unique strings
  function fms_find_unique(my_pointer, p_size) bind(c)&
  result(ntimes)
    use iso_c_binding

    type(c_ptr),            intent(in) :: my_pointer(*)     !< Array of sorted c pointer
    integer(kind=c_int),    intent(in) :: p_size            !< Size of the array
    integer(kind=c_int) :: ntimes

  end function fms_find_unique

  end interface

  !> @addtogroup fms_string_utils_mod
  !> @{
  contains

  !> @brief Converts a character array to an array of c pointers!
  !! @return An array of c pointers
  function fms_array_to_pointer(my_array) &
  result(my_pointer)
    character(len=*), target :: my_array(:) !!< Array of strings to convert
    type(c_ptr), allocatable :: my_pointer(:)

    integer :: i !< For do loops

    if (allocated(my_pointer)) call mpp_error(FATAL, "The c pointer array is &
      already allocated. Deallocated before calling fms_array_to_pointer")
    allocate(my_pointer(size(my_array)))

    do i = 1, size(my_array)
      my_pointer(i) = c_loc(my_array(i))
    enddo
  end function fms_array_to_pointer

  !> @brief Convert an array of c pointers back to a character array
  !! @return A character array
  function fms_pointer_to_array(my_pointer, narray) &
  result(my_array)
    type(c_ptr), intent(in)       :: my_pointer(*) !< Array of c pointer
    integer,     intent(in)       :: narray        !< Length of the array
    character(len=:), allocatable :: my_array(:)

    character(len=:), allocatable :: buffer !< Buffer to store a string
    integer                       :: i      !< For do loops

    allocate(character(len=255) :: my_array(narray))
    do i = 1, narray
      buffer = fms_c2f_string(my_pointer(i))
      my_array(i) = buffer
      deallocate(buffer)
    enddo
  end function fms_pointer_to_array

  !> @brief Searches through a SORTED array of pointers for a string
  !! @return the indices where the array was found
  !! If the string was not found, indices will be indices(1) = -999
  !> <br>Example usage:
  !!     my_pointer = fms_array_to_pointer(my_array)
  !!     call fms_sort_this(my_pointer, n_array, indices)
  !!     ifind = fms_find_my_string(my_pointer, n_array, string_to_find)
  function fms_find_my_string(my_pointer, narray, string_to_find) &
  result(ifind)
    type(c_ptr),      intent(in)    :: my_pointer(*)  !< Array of c pointer
    integer,          intent(in)    :: narray         !< Length of the array
    character(len=*), intent(in)    :: string_to_find !< string to find
    integer, allocatable            :: ifind(:)

    integer                       :: nfind  !< number of times the string was found
    character(len=:), allocatable :: buffer !< buffer to read the indices into

    buffer = fms_c2f_string(&
      fms_find_my_string_binding(my_pointer, narray, trim(string_to_find)//c_null_char, nfind))

    if (allocated(ifind)) call mpp_error(FATAL, "The indices array is already allocated. &
    Deallocate it before calling fms_find_my_string")

    if (nfind .gt. 0) then
      allocate(ifind(nfind))
      read(buffer,*) ifind
    else
      allocate(ifind(1))
      ifind = -999
    endif

  end function fms_find_my_string

end module fms_string_utils_mod
!> @}
! close documentation grouping
