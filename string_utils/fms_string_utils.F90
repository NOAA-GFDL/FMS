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
  use mpp_mod

  implicit none
  private

  public :: fms_array_to_pointer
  public :: fms_pointer_to_array
  public :: fms_sort_this
  public :: fms_find_my_string
  public :: fms_find_unique
  public :: fms_c2f_string
  public :: fms_f2c_string
  public :: fms_cstring2cpointer
  public :: string
  public :: string_copy
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

  !> @brief c function that finds the number of unique strings in an array of c pointers
  !! @return number of unique strings
  function fms_find_unique(my_pointer, p_size) bind(c)&
  result(ntimes)
    use iso_c_binding

    type(c_ptr),            intent(in) :: my_pointer(*)     !< Array of sorted c pointer
    integer(kind=c_int),    intent(in) :: p_size            !< Size of the array
    integer(kind=c_int) :: ntimes

  end function fms_find_unique

  !> @brief converts a kind=c_char to type c_ptr
  pure function fms_cstring2cpointer (cs) result (cp) bind(c, name="cstring2cpointer")
   import c_char, c_ptr
   character(kind=c_char), intent(in) :: cs(*) !< C string input
   type (c_ptr) :: cp !< C pointer
  end function fms_cstring2cpointer

  !> @brief Finds the length of a C-string
  integer(c_size_t) pure function c_strlen(s) bind(c,name="strlen")
    import c_size_t, c_ptr
    type(c_ptr), intent(in), value :: s !< A C-string whose size is desired
  end function

  !> @brief Frees a C pointer
  subroutine c_free(ptr) bind(c,name="free")
    import c_ptr
    type(c_ptr), value :: ptr !< A C-pointer to free
  end subroutine

end interface

!> Converts a C string to a Fortran string
!> @ingroup fms_mod
interface fms_c2f_string
  module procedure cstring_fortran_conversion
  module procedure cpointer_fortran_conversion
end interface

!> Converts a number to a string
!> @ingroup fms_mod
interface string
   module procedure string_from_integer
   module procedure string_from_real
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

  !> \brief Converts a C-string to a pointer and then to a Fortran string
  function cstring_fortran_conversion (cstring) result(fstring)
    character (kind=c_char), intent(in) :: cstring (*) !< Input C-string
    character(len=:), allocatable :: fstring    !< The fortran string returned
    fstring = cpointer_fortran_conversion(fms_cstring2cpointer(cstring))
  end function cstring_fortran_conversion

  !> \brief Converts a C-string returned from a TYPE(C_PTR) function to
  !! a fortran string with type character.
  function cpointer_fortran_conversion (cstring) result(fstring)
    type (c_ptr), intent(in) :: cstring !< Input C-pointer
    character(len=:), allocatable :: fstring    !< The fortran string returned
    character(len=:,kind=c_char), pointer :: string_buffer !< A temporary pointer to between C and Fortran
    integer(c_size_t) :: length !< The string length

    length = c_strlen(cstring)
    allocate (character(len=length, kind=c_char) :: string_buffer)
    block
      character(len=length,kind=c_char), pointer :: s
      call c_f_pointer(cstring,s)  ! Recovers a view of the C string
      string_buffer = s                   ! Copies the string contents
    end block

    allocate(character(len=length) :: fstring) !> Set the length of fstring
    fstring = string_buffer
    deallocate(string_buffer)
  end function cpointer_fortran_conversion

!> @brief Copies a Fortran string into a C string and puts c_null_char in any trailing spaces
  subroutine fms_f2c_string (dest, str_in)
    character (c_char), intent (out) :: dest (:) !< C String to be copied into
    character (len=*), intent (in) :: str_in !< Fortran string to copy to C string
    integer :: i !< for looping
!> Drop an error if the C string is not large enough to hold the input and the c_null_char at the end.
    if (len(trim(str_in)) .ge. size(dest)) call mpp_error(FATAL, &
      "The string "//trim(str_in)//" is larger than the destination C string")
!> Copy c_null_char into each spot in dest
    dest = c_null_char
!> Loop though and put each character of the Fortran string into the C string array
    do i = 1, len(trim(str_in))
      dest(i) = str_in(i:i)
    enddo
end subroutine fms_f2c_string


  !> @brief Converts an integer to a string
  !> @return The integer as a string
  function string_from_integer(i) result (res)
    integer, intent(in) :: i !< Integer to be converted to a string
    character(:),allocatable :: res !< String converted frominteger
    character(range(i)+2) :: tmp !< Temp string that is set to correct size
    write(tmp,'(i0)') i
    res = trim(tmp)
   return

  end function string_from_integer

  !#######################################################################
  !> @brief Converts a real to a string
  !> @return The real number as a string
  function string_from_real(r)
    real, intent(in) :: r !< Real number to be converted to a string
    character(len=32) :: string_from_real

    write(string_from_real,*) r

    return

  end function string_from_real

  !> @brief Safely copy a string from one buffer to another.
  subroutine string_copy(dest, source, check_for_null)
    character(len=*), intent(inout) :: dest !< Destination string.
    character(len=*), intent(in) :: source !< Source string.
    logical, intent(in), optional :: check_for_null !<Flag indicating to test for null character

    integer :: i
    logical :: check_null

    check_null = .false.
    if (present(check_for_null)) check_null = check_for_null

    i = 0
    if (check_null) then
      i = index(source, char(0)) - 1
    endif

    if (i < 1 ) i = len_trim(source)

    if (len_trim(source(1:i)) .gt. len(dest)) then
      call mpp_error(FATAL, "The input destination string is not big enough to" &
                 //" to hold the input source string.")
    endif
    dest = ""
    dest = adjustl(trim(source(1:i)))
  end subroutine string_copy

end module fms_string_utils_mod
!> @}
! close documentation grouping
