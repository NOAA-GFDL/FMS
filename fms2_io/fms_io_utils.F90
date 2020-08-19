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

!> @file

!> @brief Utility routines.
module fms_io_utils_mod
use, intrinsic :: iso_fortran_env, only: error_unit, int32, int64, real32, real64
#ifdef _OPENMP
use omp_lib
#endif
use mpp_mod
implicit none
private

public :: char_linked_list
public :: error
public :: file_exists
public :: openmp_thread_trap
public :: string_copy
public :: is_in_list
public :: append_to_list
public :: destroy_list
public :: domain_tile_filepath_mangle
public :: io_domain_tile_filepath_mangle
public :: allocate_array
public :: put_array_section
public :: get_array_section
public :: get_data_type_string
public :: get_checksum
public :: open_check
public :: string_compare
public :: restart_filepath_mangle

!> @brief A linked list of strings
type :: char_linked_list
  character(len=128) :: string
  type(char_linked_list), pointer :: head => null()
endtype char_linked_list


interface allocate_array
  module procedure allocate_array_int32_1d
  module procedure allocate_array_int32_2d
  module procedure allocate_array_int32_3d
  module procedure allocate_array_int32_4d
  module procedure allocate_array_int32_5d
  module procedure allocate_array_int64_1d
  module procedure allocate_array_int64_2d
  module procedure allocate_array_int64_3d
  module procedure allocate_array_int64_4d
  module procedure allocate_array_int64_5d
  module procedure allocate_array_real32_1d
  module procedure allocate_array_real32_2d
  module procedure allocate_array_real32_3d
  module procedure allocate_array_real32_4d
  module procedure allocate_array_real32_5d
  module procedure allocate_array_real64_1d
  module procedure allocate_array_real64_2d
  module procedure allocate_array_real64_3d
  module procedure allocate_array_real64_4d
  module procedure allocate_array_real64_5d
  module procedure allocate_array_char_1d
  module procedure allocate_array_char_2d
  module procedure allocate_array_char_3d
  module procedure allocate_array_char_4d
  module procedure allocate_array_char_5d
  module procedure allocate_array_char_6d
end interface allocate_array


interface put_array_section
  module procedure put_array_section_int32_1d
  module procedure put_array_section_int32_2d
  module procedure put_array_section_int32_3d
  module procedure put_array_section_int32_4d
  module procedure put_array_section_int32_5d
  module procedure put_array_section_int64_1d
  module procedure put_array_section_int64_2d
  module procedure put_array_section_int64_3d
  module procedure put_array_section_int64_4d
  module procedure put_array_section_int64_5d
  module procedure put_array_section_real32_1d
  module procedure put_array_section_real32_2d
  module procedure put_array_section_real32_3d
  module procedure put_array_section_real32_4d
  module procedure put_array_section_real32_5d
  module procedure put_array_section_real64_1d
  module procedure put_array_section_real64_2d
  module procedure put_array_section_real64_3d
  module procedure put_array_section_real64_4d
  module procedure put_array_section_real64_5d
end interface put_array_section


interface get_array_section
  module procedure get_array_section_int32_1d
  module procedure get_array_section_int32_2d
  module procedure get_array_section_int32_3d
  module procedure get_array_section_int32_4d
  module procedure get_array_section_int32_5d
  module procedure get_array_section_int64_1d
  module procedure get_array_section_int64_2d
  module procedure get_array_section_int64_3d
  module procedure get_array_section_int64_4d
  module procedure get_array_section_int64_5d
  module procedure get_array_section_real32_1d
  module procedure get_array_section_real32_2d
  module procedure get_array_section_real32_3d
  module procedure get_array_section_real32_4d
  module procedure get_array_section_real32_5d
  module procedure get_array_section_real64_1d
  module procedure get_array_section_real64_2d
  module procedure get_array_section_real64_3d
  module procedure get_array_section_real64_4d
  module procedure get_array_section_real64_5d
end interface get_array_section


interface get_data_type_string
  module procedure get_data_type_string_0d
  module procedure get_data_type_string_1d
  module procedure get_data_type_string_2d
  module procedure get_data_type_string_3d
  module procedure get_data_type_string_4d
  module procedure get_data_type_string_5d
end interface get_data_type_string


interface get_checksum
  module procedure get_checksum_0d
  module procedure get_checksum_1d
  module procedure get_checksum_2d
  module procedure get_checksum_3d
  module procedure get_checksum_4d
  module procedure get_checksum_5d
end interface get_checksum

contains


!> @brief Print a message to stderr, then stop the program.
subroutine error(mesg)

  character(len=*), intent(in) :: mesg !< Message that will be printed to
                                       !! stderr.

  call mpp_error(fatal, trim(mesg))
end subroutine error


!> @brief Determine if a file exists.
!! @return Flag telling if the file exists.
function file_exists(path) &
  result(exists)

  character(len=*), intent(in) :: path !< Path to file.
  logical :: exists

!$omp critical (file_existence_inquiry)
  inquire(file=trim(path), exist=exists)
!$omp end critical (file_existence_inquiry)
end function file_exists


!> @brief Catch OpenMP parallel machines.
subroutine openmp_thread_trap()

#ifdef _OPENMP
    if (omp_get_level() .gt. 0) then
      call error("this routine is not thread-safe.  Please do not" &
                 //" call it in an OpenMP threaded region.")
    endif
#endif
end subroutine openmp_thread_trap


!> @brief Safely copy a string from one buffer to another.
subroutine string_copy(dest, source, check_for_null)
  character(len=*), intent(inout) :: dest !< Destination string.
  character(len=*), intent(in) :: source !< Source string.
  logical, intent(in), optional :: check_for_null !<Flag indicating to test for null character

  integer :: i
  logical :: check_null

  check_null = .false.
  if (present(check_for_null)) check_null = check_for_null
  if (len_trim(source) .gt. len(dest)) then
    call error("The input destination string is not big enough to" &
                 //" to hold the input source string.")
  endif
  dest = ""
  dest = adjustl(trim(source))

  if (check_null) then
     i = 0
     i = index(dest, char(0))
     if (i > 0 ) dest = dest(1:i-1)
  endif

end subroutine string_copy


!> @brief Compare strings.
!! @return Flag telling if the strings are the same.
function string_compare(string1, string2, ignore_case) &
  result(same)

  character(len=*), intent(in) :: string1
  character(len=*), intent(in) :: string2
  logical,intent(in), optional :: ignore_case
  logical :: same

  if (len_trim(string1) .ne. len_trim(string2)) then
    same = .false.
    return
  endif
  if (present(ignore_case)) then
    if (ignore_case) then
      same = trim(lowercase(string1)) .eq. trim(lowercase(string2))
      return
    endif
  endif
  same = trim(string1) .eq. trim(string2)
end function string_compare


!> @brief Determine if a string exists in a character linked list.
!! @return Flag telling if the string is found in the list.
function is_in_list(list, string, ignore_case) &
  result(in_list)

  type(char_linked_list), pointer, intent(in) :: list !< Linked list.
  character(len=*), intent(in) :: string !< Input string.
  logical, intent(in), optional :: ignore_case !< Flag to ignore case.
  logical :: in_list

  type(char_linked_list), pointer :: p

  in_list = .false.
  p => list
  do while(associated(p))
    if (string_compare(string, p%string, ignore_case)) then
      in_list = .true.
      return
    else
      p => p%head
    endif
  enddo
end function is_in_list


!> @brief Add node to character linked list.
subroutine append_to_list(list, string)
  type(char_linked_list), pointer, intent(inout) :: list !< Linked list.
  character(len=*), intent(in) :: string !< Input string.

  type(char_linked_list), pointer :: node
  type(char_linked_list), pointer :: p

  allocate(node)
  call string_copy(node%string, string)
  node%head => null()
  if (associated(list)) then
    p => list
    do while (associated(p%head))
      p => p%head
    enddo
    p%head => node
  else
    list => node
  endif
end subroutine append_to_list


!> @brief Deallocate all nodes on a character linked list.
subroutine destroy_list(list)

  type(char_linked_list), pointer, intent(inout) :: list !< Linked list.

  type(char_linked_list), pointer :: p
  type(char_linked_list), pointer :: p2
  p => list
  do while (associated(p))
    p2 => p%head
    deallocate(p)
    p => p2
  enddo
  list => null()
end subroutine destroy_list


!> @brief Determine if the "domain tile string" (.tilex.) exists in the input filename.
!! @internal
function has_domain_tile_string(string) &
  result(has_string)

  character(len=*), intent(in) :: string !< Input string.
  logical :: has_string

  integer :: l
  integer :: i, j

  has_string = .false.
! Assigns i to the index where ".tile" starts
  i = index(trim(string), ".tile", back=.true.)
  if (i .ne. 0) then
    l = len_trim(string)
! Sets i to the index after .tile
    i = i + 5
    j = i
    do while (i .le. l)
! If the ith characters is a dot but i not equal to the index after .tile set has_string to true
      if (verify(string(i:i), ".") .eq. 0 .and. j .ne. i) then
        has_string = .true.
        exit
! If the ith characters is NOT a number exit function and has_string will stay as false
      elseif (verify(string(i:i), "0123456789") .ne. 0) then
        exit
      endif
      i = i + 1
    enddo
  endif
end function has_domain_tile_string


!> @brief Add the domain tile id to an input filepath.
!! @internal
subroutine domain_tile_filepath_mangle(dest, source, domain_tile_id)

  character(len=*), intent(inout) :: dest !< Output filepath.
  character(len=*), intent(in) :: source !< Input filepath.
  integer, intent(in) :: domain_tile_id !< Domain tile id.

  integer :: i

  if (has_domain_tile_string(source)) then
    call error("this file has already had a domain tile id added.")
  endif
  i = index(trim(source), ".nc", back=.true.)
  if (i .eq. 0) then
    call error("file "//trim(source)//" does not contain .nc")
  endif
  write(dest, '(a,i1,a)') source(1:i-1)//".tile", &
                          domain_tile_id, source(i:len_trim(source))
end subroutine domain_tile_filepath_mangle


!> @brief Determine if the "I/O domain tile string" (.nc.xxxx) exists in the input filename.
!! @internal
function has_io_domain_tile_string(string) &
  result (has_string)

  character(len=*), intent(in) :: string !< Input string.
  logical :: has_string

  integer :: i
  integer :: l

  has_string = .false.
  i = index(trim(string), ".nc.", back=.true.)
  if (i .ne. 0) then
    l = len_trim(string)
    i = i + 1
    do while (i .le. l)
      if (verify(string(i:i), "0123456789") .ne. 0) then
        return
      endif
      i = i + 1
    enddo
    has_string = .true.
  endif
end function has_io_domain_tile_string


!> @brief Add the I/O domain tile id to an input filepath.
!!
subroutine io_domain_tile_filepath_mangle(dest, source, io_domain_tile_id)

  character(len=*), intent(inout) :: dest !< Output filepath.
  character(len=*), intent(in) :: source !< Input filepath.
  integer, intent(in) :: io_domain_tile_id !< I/O domain tile id.

  if (has_io_domain_tile_string(source)) then
    call error("this file has already had a domain tile id added.")
  endif
  write(dest,'(a,i4.4)') trim(source)//".", io_domain_tile_id
end subroutine io_domain_tile_filepath_mangle


!> @brief Determine if the "restart string" (.res.) exists in the input filename.
!! @internal
function has_restart_string(string) &
  result (has_string)

  character(len=*), intent(in) :: string !< Input string.
  logical :: has_string

  has_string = index(trim(string), ".res.", back=.true.) .ne. 0
end function has_restart_string


!> @brief Add ".res" to an input file path.
!!
subroutine restart_filepath_mangle(dest, source)

  character(len=*), intent(inout) :: dest
  character(len=*), intent(in) :: source

  integer :: i

  if (has_restart_string(source)) then
    call string_copy(dest, source)
    return
  endif
  if (has_domain_tile_string(source)) then
    i = index(trim(source), ".tile", back=.true.)
  else
    i = index(trim(source), ".nc", back=.true.)
    if (i .eq. 0) then
      call error("file "//trim(source)//" does not contain .nc")
    endif
  endif
  call string_copy(dest, source(1:i-1)//".res"//source(i:len_trim(source)))
end subroutine restart_filepath_mangle

subroutine open_check(flag, fname)

  logical, intent(in) :: flag
  character(len=*), intent(in), optional :: fname !< The file name

  if (.not. flag) then
     if (present(fname)) then
          call mpp_error(fatal, "Error occured while opening file "//trim(fname))
     else
          call mpp_error(fatal, "Error occured while opening file.")
     endif
  endif
end subroutine open_check


include "array_utils.inc"
include "array_utils_char.inc"
include "get_data_type_string.inc"
include "get_checksum.inc"


end module fms_io_utils_mod
