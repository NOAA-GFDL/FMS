!> @file

!> @brief Utility routines.
module fms_io_utils_mod
use, intrinsic :: iso_fortran_env, only: error_unit, int32, int64, real32, real64
#ifdef _OPENMP
use omp_lib
#endif
use mpp_mod


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


!> @brief A linked list of strings.
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
subroutine string_copy(dest, source)
  character(len=*), intent(inout) :: dest !< Destination string.
  character(len=*), intent(in) :: source !< Source string.

  if (len_trim(source) .gt. len(dest)) then 
    call error("The input destination string is not big enough to" &
                 //" to hold the input source string.")
  endif
  dest = ""
  dest = adjustl(trim(source))
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


!> @brief Add the domain tile id to an input filepath.
!! @internal
subroutine domain_tile_filepath_mangle(dest, source, domain_tile_id)
  character(len=*), intent(inout) :: dest !< Output filepath.
  character(len=*), intent(in) :: source !< Input filepath.
  integer, intent(in) :: domain_tile_id !< Domain tile id.

  integer :: i

  i = index(trim(source), ".nc", back=.true.)
  if (i .eq. 0) then
    call error("file "//trim(source)//" does not contain .nc")
  endif
  write(dest, '(a,i1,a)') source(1:i-1)//".tile", &
                          domain_tile_id, source(i:len_trim(source))
end subroutine domain_tile_filepath_mangle


!> @brief Add the I/O domain tile id to an input filepath.
!! @internal
subroutine io_domain_tile_filepath_mangle(dest, source, io_domain_tile_id)

  character(len=*), intent(inout) :: dest !< Output filepath.
  character(len=*), intent(in) :: source !< Input filepath.
  integer, intent(in) :: io_domain_tile_id !< I/O domain tile id.

  write(dest,'(a,i4.4)') trim(source)//".", io_domain_tile_id
end subroutine io_domain_tile_filepath_mangle


!> @brief Add ".res" to an input file path.
!! @internal
subroutine restart_filepath_mangle(dest, source)

  character(len=*), intent(inout) :: dest
  character(len=*), intent(in) :: source

  character(len=512) :: buf
  integer :: i
  integer :: j
  integer :: k

  if (index(trim(source), ".res.", back=.true.) .ne. 0) then
    call string_copy(dest, source)
    return
  endif
  i = index(trim(source), ".tile", back=.true.)
  if (i .eq. 0) then
    i = index(trim(source), ".nc", back=.true.)
    if (i .eq. 0) then
      call error("file "//trim(source)//" does not contain .nc")
    endif
  else
    buf = trim(source(i+5:len(source)))
    j = index(trim(buf),".")
    if (j .eq. 0) then
      call error("file "//trim(source)//" does not contain .tilex.")
    endif
    do k = 1, j-1
      if (verify(buf(k:k), "0123456789") .ne. 0) then
        call error("file "//trim(source)//" does not contain .tilex.")
      endif
    enddo
  endif
  call string_copy(dest, source(1:i-1)//".res"//source(i:len_trim(source)))
end subroutine restart_filepath_mangle


!> @brief Create a new file path.
!! @internal
subroutine get_new_filename(path, new_path, directory, timestamp, new_name)

  character(len=*), intent(in) :: path !< File path.
  character(len=*), intent(out) :: new_path !< New file path.
  character(len=*), intent(in), optional :: directory !< Directory
  character(len=*), intent(in), optional :: timestamp !< Time.
  character(len=*), intent(in), optional :: new_name !< New file basename.

  character(len=256) :: dir
  character(len=256) :: tstamp
  character(len=256) :: nname

  dir = ""
  if (present(directory)) then
    call string_copy(dir, trim(directory)//"/")
  endif
  tstamp = ""
  if (present(timestamp)) then
    call string_copy(tstamp, timestamp//".")
  endif
  call string_copy(nname, trim(path))
  if (present(new_name)) then
    call string_copy(nname, new_name)
  endif
  call string_copy(new_path, trim(dir)//trim(tstamp)//trim(nname))
end subroutine get_new_filename


include "array_utils.inc"
include "array_utils_char.inc"
include "get_data_type_string.inc"
include "get_checksum.inc"


end module fms_io_utils_mod
