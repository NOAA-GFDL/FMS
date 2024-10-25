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
!> @defgroup fm_util_mod fm_util_mod
!> @ingroup field_manager
!> @brief This module provides utility routines for the field manager.
!!
!#######################################################################

!> Set a real array in the Field Manager tree.
subroutine FM_UTIL_SET_VALUE_REAL_ARRAY_(name, rval, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
real(FMS_FM_KIND_), intent(in)                          :: rval(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=FMS_PATH_LEN), intent(in), optional       :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_real_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=FMS_PATH_LEN)     :: good_name_list_use
logical                         :: add_name

integer, parameter :: lkind=FMS_FM_KIND_

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, 0.0_lkind, index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, rval(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, rval(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, rval(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine FM_UTIL_SET_VALUE_REAL_ARRAY_  !}

!#######################################################################

!> Set a real value in the Field Manager tree.
subroutine FM_UTIL_SET_VALUE_REAL_(name, rval, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
real(FMS_FM_KIND_), intent(in)          :: rval
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_real'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=FMS_PATH_LEN)     :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // &
                   &  ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, rval, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, rval, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, rval, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, rval)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, rval)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine FM_UTIL_SET_VALUE_REAL_  !}

!#######################################################################
!> @}
! close documentation grouping
