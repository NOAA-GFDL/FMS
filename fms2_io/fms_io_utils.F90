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
!> @defgroup fms_io_utils_mod fms_io_utils_mod
!> @ingroup fms2_io
!> @brief Misc. utility routines for use in @ref fms2_io

!> @file
!> @brief File for @ref fms_io_utils_mod

!> @addtogroup fms_io_utils_mod
!> @{
module fms_io_utils_mod
use, intrinsic :: iso_fortran_env, only: error_unit
!use mpp_mod, only : get_ascii_file_num_lines_and_length, read_ascii_file
#ifdef _OPENMP
use omp_lib
#endif
use mpp_mod
use mpp_domains_mod, only: domain2D, domainUG, mpp_get_ntile_count, &
                           mpp_get_current_ntile, mpp_get_tile_id, &
                           mpp_get_UG_domain_ntiles, mpp_get_UG_domain_tile_id
use platform_mod
implicit none
private

character(len=32), save :: filename_appendix = '' !< Appendix added to the restart filename

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
public :: string2
public :: open_check
public :: string_compare
public :: restart_filepath_mangle
public :: ascii_read
public :: parse_mask_table
public :: get_mosaic_tile_file
public :: get_filename_appendix
public :: set_filename_appendix
public :: get_instance_filename
public :: nullify_filename_appendix

!> @}

!> @brief A linked list of strings
!> @ingroup fms_io_utils_mod
type :: char_linked_list
  character(len=128) :: string
  type(char_linked_list), pointer :: head => null()
endtype char_linked_list

!> @brief Converts a given integer or real into a character string
!> @ingroup fms_io_utils_mod
interface string2
  module procedure string_from_integer2
  module procedure string_from_real2
end interface string2


!> @ingroup fms_io_utils_mod
interface parse_mask_table
  module procedure parse_mask_table_2d
  module procedure parse_mask_table_3d
end interface parse_mask_table

!> @ingroup fms_io_utils_mod
interface get_mosaic_tile_file
  module procedure get_mosaic_tile_file_sg
  module procedure get_mosaic_tile_file_ug
end interface get_mosaic_tile_file

!> @ingroup fms_io_utils_mod
interface allocate_array
  module procedure allocate_array_i4_kind_1d
  module procedure allocate_array_i4_kind_2d
  module procedure allocate_array_i4_kind_3d
  module procedure allocate_array_i4_kind_4d
  module procedure allocate_array_i4_kind_5d
  module procedure allocate_array_i8_kind_1d
  module procedure allocate_array_i8_kind_2d
  module procedure allocate_array_i8_kind_3d
  module procedure allocate_array_i8_kind_4d
  module procedure allocate_array_i8_kind_5d
  module procedure allocate_array_r4_kind_1d
  module procedure allocate_array_r4_kind_2d
  module procedure allocate_array_r4_kind_3d
  module procedure allocate_array_r4_kind_4d
  module procedure allocate_array_r4_kind_5d
  module procedure allocate_array_r8_kind_1d
  module procedure allocate_array_r8_kind_2d
  module procedure allocate_array_r8_kind_3d
  module procedure allocate_array_r8_kind_4d
  module procedure allocate_array_r8_kind_5d
  module procedure allocate_array_char_1d
  module procedure allocate_array_char_2d
  module procedure allocate_array_char_3d
  module procedure allocate_array_char_4d
  module procedure allocate_array_char_5d
  module procedure allocate_array_char_6d
end interface allocate_array


!> @ingroup fms_io_utils_mod
interface put_array_section
  module procedure put_array_section_i4_kind_1d
  module procedure put_array_section_i4_kind_2d
  module procedure put_array_section_i4_kind_3d
  module procedure put_array_section_i4_kind_4d
  module procedure put_array_section_i4_kind_5d
  module procedure put_array_section_i8_kind_1d
  module procedure put_array_section_i8_kind_2d
  module procedure put_array_section_i8_kind_3d
  module procedure put_array_section_i8_kind_4d
  module procedure put_array_section_i8_kind_5d
  module procedure put_array_section_r4_kind_1d
  module procedure put_array_section_r4_kind_2d
  module procedure put_array_section_r4_kind_3d
  module procedure put_array_section_r4_kind_4d
  module procedure put_array_section_r4_kind_5d
  module procedure put_array_section_r8_kind_1d
  module procedure put_array_section_r8_kind_2d
  module procedure put_array_section_r8_kind_3d
  module procedure put_array_section_r8_kind_4d
  module procedure put_array_section_r8_kind_5d
end interface put_array_section


!> @ingroup fms_io_utils_mod
interface get_array_section
  module procedure get_array_section_i4_kind_1d
  module procedure get_array_section_i4_kind_2d
  module procedure get_array_section_i4_kind_3d
  module procedure get_array_section_i4_kind_4d
  module procedure get_array_section_i4_kind_5d
  module procedure get_array_section_i8_kind_1d
  module procedure get_array_section_i8_kind_2d
  module procedure get_array_section_i8_kind_3d
  module procedure get_array_section_i8_kind_4d
  module procedure get_array_section_i8_kind_5d
  module procedure get_array_section_r4_kind_1d
  module procedure get_array_section_r4_kind_2d
  module procedure get_array_section_r4_kind_3d
  module procedure get_array_section_r4_kind_4d
  module procedure get_array_section_r4_kind_5d
  module procedure get_array_section_r8_kind_1d
  module procedure get_array_section_r8_kind_2d
  module procedure get_array_section_r8_kind_3d
  module procedure get_array_section_r8_kind_4d
  module procedure get_array_section_r8_kind_5d
end interface get_array_section


!> @ingroup fms_io_utils_mod
interface get_data_type_string
  module procedure get_data_type_string_0d
  module procedure get_data_type_string_1d
  module procedure get_data_type_string_2d
  module procedure get_data_type_string_3d
  module procedure get_data_type_string_4d
  module procedure get_data_type_string_5d
end interface get_data_type_string

!> @addtogroup fms_io_utils_mod
!> @{
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

  i = 0
  if (check_null) then
     i = index(source, char(0)) - 1
  endif

  if (i < 1 ) i = len_trim(source)

  if (len_trim(source(1:i)) .gt. len(dest)) then
    call error("The input destination string is not big enough to" &
                 //" to hold the input source string.")
  endif
  dest = ""
  dest = adjustl(trim(source(1:i)))

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
    call error("The file "//trim(source)//" has a domain tile id (tileX) added. Check your open_file call")
  endif
  i = index(trim(source), ".nc", back=.true.)
  if (i .eq. 0) then
    call error("The file "//trim(source)//" does not contain .nc. Check your open_file call")
  endif
  write(dest, '(a,i0,a)') source(1:i-1)//".tile", &
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
subroutine io_domain_tile_filepath_mangle(dest, source, io_domain_tile_id)

  character(len=*), intent(inout) :: dest !< Output filepath.
  character(len=*), intent(in) :: source !< Input filepath.
  integer, intent(in) :: io_domain_tile_id !< I/O domain tile id.

  if (has_io_domain_tile_string(source)) then
    call error("The file "//trim(source)//" has already had a domain tile id (.nc.XXXX) added. Check your open_file call.")
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
      call error("The file "//trim(source)//" does not contain .nc. Check your open_file call")
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

!> @brief Read the ascii text from filename `ascii_filename`into string array
!! `ascii_var`
subroutine ascii_read(ascii_filename, ascii_var, num_lines, max_length)
  character(len=*), intent(in) :: ascii_filename !< The file name to be read
  character(len=:), dimension(:), allocatable, intent(out) :: ascii_var !< The
                                                                        !! string
                                                                        !! array
  integer, optional, intent(out) :: num_lines !< Optional argument to return number of lines in file
  integer, optional, intent(out) :: max_length !< Optional argument to return max_length of line in file
  integer, dimension(2) :: lines_and_length !< lines = 1, length = 2
  if(allocated(ascii_var)) deallocate(ascii_var)
  lines_and_length = get_ascii_file_num_lines_and_length(ascii_filename)
  allocate(character(len=lines_and_length(2))::ascii_var(lines_and_length(1)))
  call read_ascii_file(ascii_filename, lines_and_length(2), ascii_var)
  if(present(num_lines)) num_lines = lines_and_length(1)
  if(present(max_length)) max_length = lines_and_length(2)
end subroutine ascii_read

!> @brief Populate 2D maskmap from mask_table given a model
subroutine parse_mask_table_2d(mask_table, maskmap, modelname)

  character(len=*), intent(in) :: mask_table !< Mask table to be read in
  logical,         intent(out) :: maskmap(:,:) !< 2D Mask output
  character(len=*), intent(in) :: modelname !< Model to which this applies

  integer                      :: nmask, layout(2)
  integer, allocatable         :: mask_list(:,:)
  character(len=:), dimension(:), allocatable :: mask_table_contents
  integer                      :: iocheck, n, stdoutunit, offset
  character(len=128)           :: record

  maskmap = .true.
  nmask = 0
  stdoutunit = stdout()
  call ascii_read(mask_table, mask_table_contents)
  if( mpp_pe() == mpp_root_pe() ) then
     read(mask_table_contents(1), FMT=*, IOSTAT=iocheck) nmask
     if (iocheck > 0) then
         call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error in reading nmask from file variable")
     elseif (iocheck < 0) then
         call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error: nmask not completely read from file variable")
     endif
     write(stdoutunit,*)"parse_mask_table: Number of domain regions masked in ", trim(modelname), " = ", nmask
     if( nmask > 0 ) then
        !--- read layout from mask_table and confirm it matches the shape of maskmap
        read(mask_table_contents(2), FMT=*, IOSTAT=iocheck) layout
        if (iocheck > 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error in reading layout from file variable")
        elseif (iocheck < 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error: layout not completely read from file variable")
        endif
        if( (layout(1) .NE. size(maskmap,1)) .OR. (layout(2) .NE. size(maskmap,2)) )then
           write(stdoutunit,*)"layout=", layout, ", size(maskmap) = ", size(maskmap,1), size(maskmap,2)
           call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): layout in file "//trim(mask_table)// &
                  "does not match size of maskmap for "//trim(modelname))
        endif
        !--- make sure mpp_npes() == layout(1)*layout(2) - nmask
        if( mpp_npes() .NE. layout(1)*layout(2) - nmask ) call mpp_error(FATAL, &
           "fms2_io(parse_mask_table_2d): mpp_npes() .NE. layout(1)*layout(2) - nmask for "//trim(modelname))
     endif
   endif

   call mpp_broadcast(nmask, mpp_root_pe())

   if(nmask==0) return

   allocate(mask_list(nmask,2))

   if( mpp_pe() == mpp_root_pe() ) then
     n = 0
     offset = 3
     do while (offset + n < size(mask_table_contents)+1)
        read(mask_table_contents(n+offset),'(a)',iostat=iocheck) record
        if (iocheck > 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error in reading record from file variable")
        elseif (iocheck < 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error: record not completely read from file variable")
        endif
        if (record(1:1) == '#') then
            offset = offset + 1
            cycle
        elseif (record(1:10) == '          ') then
            offset = offset + 1
            cycle
        endif
        n = n + 1
        if( n > nmask ) then
           call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): number of mask_list entry "// &
                "is greater than nmask in file "//trim(mask_table) )
        endif
        read(record,*,iostat=iocheck) mask_list(n,1), mask_list(n,2)
        if (iocheck > 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error in reading mask_list from record variable")
        elseif (iocheck < 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_2d): Error: mask_list not completely read from record variable")
        endif
     enddo

     !--- make sure the number of entry for mask_list is nmask
     if( n .NE. nmask) call mpp_error(FATAL, &
        "fms2_io(parse_mask_table_2d): number of mask_list entry does not match nmask in file "//trim(mask_table))
  endif

  call mpp_broadcast(mask_list, 2*nmask, mpp_root_pe())
  do n = 1, nmask
     maskmap(mask_list(n,1),mask_list(n,2)) = .false.
  enddo

  deallocate(mask_list)

end subroutine parse_mask_table_2d


!> @brief Populate 3D maskmap from mask_table given a model
subroutine parse_mask_table_3d(mask_table, maskmap, modelname)

  character(len=*), intent(in) :: mask_table !< Mask table to be read in
  logical,         intent(out) :: maskmap(:,:,:) !< 2D Mask output
  character(len=*), intent(in) :: modelname !< Model to which this applies

  integer                      :: nmask, layout(2)
  integer, allocatable         :: mask_list(:,:)
  character(len=:), dimension(:), allocatable :: mask_table_contents
  integer                      :: iocheck, n, stdoutunit, ntiles, offset
  character(len=128)           :: record

  maskmap = .true.
  nmask = 0
  stdoutunit = stdout()
  call ascii_read(mask_table, mask_table_contents)
  if( mpp_pe() == mpp_root_pe() ) then
     read(mask_table_contents(1), FMT=*, IOSTAT=iocheck) nmask
     if (iocheck > 0) then
         call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error in reading nmask from file variable")
     elseif (iocheck < 0) then
         call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error: nmask not completely read from file variable")
     endif
     write(stdoutunit,*)"parse_mask_table: Number of domain regions masked in ", trim(modelname), " = ", nmask
     if( nmask > 0 ) then
        !--- read layout from mask_table and confirm it matches the shape of maskmap
        read(mask_table_contents(2), FMT=*, IOSTAT=iocheck) layout(1), layout(2), ntiles
        if (iocheck > 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error in reading layout from file variable")
        elseif (iocheck < 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error: layout not completely read from file variable")
        endif
        if( (layout(1) .NE. size(maskmap,1)) .OR. (layout(2) .NE. size(maskmap,2)) )then
           write(stdoutunit,*)"layout=", layout, ", size(maskmap) = ", size(maskmap,1), size(maskmap,2)
           call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): layout in file "//trim(mask_table)// &
                  "does not match size of maskmap for "//trim(modelname))
        endif
        if( ntiles .NE. size(maskmap,3) ) then
           write(stdoutunit,*)"ntiles=", ntiles, ", size(maskmap,3) = ", size(maskmap,3)
           call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): ntiles in file "//trim(mask_table)// &
                  "does not match size of maskmap for "//trim(modelname))
        endif
        !--- make sure mpp_npes() == layout(1)*layout(2) - nmask
        if( mpp_npes() .NE. layout(1)*layout(2)*ntiles - nmask ) then
           print*, "layout=", layout, nmask, mpp_npes()
           call mpp_error(FATAL, &
              "fms2_io(parse_mask_table_3d): mpp_npes() .NE. layout(1)*layout(2) - nmask for "//trim(modelname))
        endif
      endif
   endif

   call mpp_broadcast(nmask, mpp_root_pe())

   if(nmask==0) return

   allocate(mask_list(nmask,3))

   if( mpp_pe() == mpp_root_pe() ) then
     n = 0
     offset = 3
     do while (offset + n < size(mask_table_contents)+1)
        read(mask_table_contents(n+offset),'(a)',iostat=iocheck) record
        if (iocheck > 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error in reading record from file variable")
        elseif (iocheck < 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error: record not completely read from file variable")
        endif
        if (record(1:1) == '#') then
            offset = offset + 1
            cycle
        elseif (record(1:10) == '          ') then
            offset = offset + 1
            cycle
        endif
        n = n + 1
        if( n > nmask ) then
           call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): number of mask_list entry "// &
                "is greater than nmask in file "//trim(mask_table) )
        endif
        read(record,*,iostat=iocheck) mask_list(n,1), mask_list(n,2), mask_list(n,3)
        if (iocheck > 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error in reading mask_list from record variable")
        elseif (iocheck < 0) then
            call mpp_error(FATAL, "fms2_io(parse_mask_table_3d): Error: mask_list not completely read from record variable")
        endif
     enddo

     !--- make sure the number of entry for mask_list is nmask
     if( n .NE. nmask) call mpp_error(FATAL, &
        "fms2_io(parse_mask_table_3d): number of mask_list entry does not match nmask in file "//trim(mask_table))
!     call mpp_close(unit)
  endif

  call mpp_broadcast(mask_list, 3*nmask, mpp_root_pe())
  do n = 1, nmask
     maskmap(mask_list(n,1),mask_list(n,2),mask_list(n,3)) = .false.
  enddo

  deallocate(mask_list)
end subroutine parse_mask_table_3d

!> @brief Determine tile_file for structured grid based on filename and current
!! tile on mpp_domain (this is mostly used for ongrid data_overrides)
subroutine get_mosaic_tile_file_sg(file_in, file_out, is_no_domain, domain, tile_count)
  character(len=*), intent(in)            :: file_in !< name of 'base' file
  character(len=*), intent(out)           :: file_out !< name of tile_file
  logical,          intent(in)            :: is_no_domain !< are we providing a
                                                          !! domain
  type(domain2D),   intent(in), optional, target :: domain !< domain provided
  integer,          intent(in), optional  :: tile_count !< tile count

  character(len=256)                             :: basefile, tilename
  character(len=2)                               :: my_tile_str
  integer                                        :: lens, ntiles, ntileMe, tile, my_tile_id
  integer, dimension(:), allocatable             :: tile_id
  type(domain2d), pointer, save                  :: d_ptr =>NULL()
  logical                                        :: domain_exist

  if(index(file_in, '.nc', back=.true.)==0) then
     basefile = trim(file_in)
  else
     lens = len_trim(file_in)
     if(file_in(lens-2:lens) .NE. '.nc') call mpp_error(FATAL, &
          'get_mosaic_tile_file_sg: .nc should be at the end of file '//trim(file_in))
     basefile = file_in(1:lens-3)
  end if

  !--- get the tile name
  ntiles = 1
  my_tile_id = 1
  domain_exist = .false.
  if(PRESENT(domain))then
     domain_exist = .true.
     ntiles = mpp_get_ntile_count(domain)
     d_ptr => domain
  endif

  if(domain_exist) then
     ntileMe = mpp_get_current_ntile(d_ptr)
     allocate(tile_id(ntileMe))
     tile_id = mpp_get_tile_id(d_ptr)
     tile = 1
     if(present(tile_count)) tile = tile_count
     my_tile_id = tile_id(tile)
  endif

  if(ntiles > 1 .or. my_tile_id > 1 )then
     write(my_tile_str, '(I0)') my_tile_id
     tilename = 'tile'//trim(my_tile_str)
     if(index(basefile,'.'//trim(tilename),back=.true.) == 0)then
        basefile = trim(basefile)//'.'//trim(tilename);
     end if
  end if
  if(allocated(tile_id)) deallocate(tile_id)

  file_out = trim(basefile)//'.nc'

  d_ptr =>NULL()

end subroutine get_mosaic_tile_file_sg

!> @brief Determine tile_file for unstructured grid based on filename and current
!! tile on mpp_domain (this is mostly used for ongrid data_overrides)
subroutine get_mosaic_tile_file_ug(file_in, file_out, domain)
  character(len=*), intent(in)           :: file_in !< name of base file
  character(len=*), intent(out)          :: file_out !< name of tile file
  type(domainUG),   intent(in), optional :: domain !< domain provided

  character(len=256)                     :: basefile, tilename
  character(len=2)                       :: my_tile_str
  integer                                :: lens, ntiles, my_tile_id

  if(index(file_in, '.nc', back=.true.)==0) then
     basefile = trim(file_in)
  else
     lens = len_trim(file_in)
     if(file_in(lens-2:lens) .NE. '.nc') call mpp_error(FATAL, &
          'fms_io_mod: .nc should be at the end of file '//trim(file_in))
     basefile = file_in(1:lens-3)
  end if

  !--- get the tile name
  ntiles = 1
  my_tile_id = 1
  if(PRESENT(domain))then
     ntiles = mpp_get_UG_domain_ntiles(domain)
     my_tile_id = mpp_get_UG_domain_tile_id(domain)
  endif

  if(ntiles > 1 .or. my_tile_id > 1 )then
     write(my_tile_str, '(I0)') my_tile_id
     tilename = 'tile'//trim(my_tile_str)
     if(index(basefile,'.'//trim(tilename),back=.true.) == 0)then
        basefile = trim(basefile)//'.'//trim(tilename);
     end if
  end if

  file_out = trim(basefile)//'.nc'

end subroutine get_mosaic_tile_file_ug

!> @brief Writes filename appendix to "string_out"
subroutine get_filename_appendix(string_out)
  character(len=*) , intent(out) :: string_out !< String to write the filename_appendix to

  string_out = trim(filename_appendix)

end subroutine get_filename_appendix

!> @brief Clears the filename_appendix module variable
subroutine nullify_filename_appendix()

  filename_appendix = ''

end subroutine nullify_filename_appendix

!> @brief Save "string_in" as a module variable that will added to the filename
!! of the restart files
subroutine set_filename_appendix(string_in)
  character(len=*) , intent(in) :: string_in !< String that will be saved as a module variable

  ! Check if string has already been added
  if (len_trim(filename_appendix) > 0) then
      call error("Set_filename_appendix: The filename appendix has already be set " &
                 //"call 'nullify_filename_appendix' first")
  endif

  filename_appendix = trim(string_in)

end subroutine set_filename_appendix

!> @brief Adds the filename_appendix to name_in and sets it as name_out
subroutine get_instance_filename(name_in,name_out)
  character(len=*)  , intent(in)  :: name_in  !< Buffer to add the filename_appendix to
  character(len=*), intent(inout) :: name_out !< name_in with the filename_appendix

  integer :: length !< Length of name_in
  integer :: i !< no description

  length = len_trim(name_in)
  name_out = name_in(1:length)

  if(len_trim(filename_appendix) > 0) then
     !< If .tileXX is in the filename add the appendix before it
     if (has_domain_tile_string(name_in)) then
         i = index(trim(name_in), ".tile", back=.true.)
         name_out = name_in(1:i-1)    //'.'//trim(filename_appendix)//name_in(i:length)
         return
     endif

     !< If .nc is in the filename add the appendix before it
     i = index(trim(name_in), ".nc", back=.true.)
     if ( i .ne. 0 ) then
        name_out = name_in(1:i-1)//'.'//trim(filename_appendix)//name_in(i:length)
     else
     !< If .nc is not in the name, add the appendix at the end of the file
        name_out = name_in(1:length)  //'.'//trim(filename_appendix)
     end if
  end if

end subroutine get_instance_filename

function string_from_integer2(n)
    integer, intent(in) :: n
    character(len=16) :: string_from_integer2
    if(n<0) then
       call mpp_error(FATAL, 'fms2_io_mod: n should be non-negative integer, contact developer')
    else if( n<10 ) then
       write(string_from_integer2,'(i1)') n
    else if( n<100 ) then
       write(string_from_integer2,'(i2)') n
    else if( n<1000 ) then
       write(string_from_integer2,'(i3)') n
    else if( n<10000 ) then
       write(string_from_integer2,'(i4)') n
    else if( n<100000 ) then
       write(string_from_integer2,'(i5)') n
    else if( n<1000000 ) then
       write(string_from_integer2,'(i6)') n
    else if( n<10000000 ) then
       write(string_from_integer2,'(i7)') n
    else if( n<100000000 ) then
       write(string_from_integer2,'(i8)') n
    else
       call mpp_error(FATAL, 'fms2_io_mod: n is greater than 1e8, contact developer')
    end if

    return

end function string_from_integer2

function string_from_real2(a)
    real, intent(in) :: a
    character(len=32) :: string_from_real2

    write(string_from_real2,*) a

    return

end function string_from_real2

include "array_utils.inc"
include "array_utils_char.inc"
include "get_data_type_string.inc"


end module fms_io_utils_mod
!> @}
! close documentation grouping
