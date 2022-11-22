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

!> @defgroup yaml_parser_mod yaml_parser_mod
!> @ingroup parser
!> @brief Routines to use for parsing yaml files

!> @file
!> @brief File for @ref yaml_parser_mod

!> @addtogroup yaml_parser_mod
!> @{
module yaml_parser_mod

#ifdef use_yaml
use fms_mod, only: fms_c2f_string
use fms_string_utils_mod, only: string_copy
use platform_mod
use mpp_mod
use iso_c_binding

implicit none
private

public :: open_and_parse_file
public :: get_num_blocks
public :: get_block_ids
public :: get_value_from_key
public :: get_nkeys
public :: get_key_ids
public :: get_key_name
public :: get_key_value
!public :: clean_up
!> @}

!> @brief Dermine the value of a key from a keyname
!> @ingroup yaml_parser_mod
interface get_value_from_key
  module procedure get_value_from_key_0d
  module procedure get_value_from_key_1d
end interface get_value_from_key

!> @brief c functions binding
!> @ingroup yaml_parser_mod
interface

!> @brief Private c function that opens and parses a yaml file (see yaml_parser_binding.c)
!! @return Flag indicating if the read was sucessful
function open_and_parse_file_wrap(filename, file_id) bind(c) &
   result(sucess)
   use iso_c_binding, only: c_char, c_int, c_bool
   character(kind=c_char), intent(in) :: filename(*) !< Filename of the yaml file
   integer(kind=c_int), intent(out) :: file_id !< File id corresponding to the yaml file that was opened
   logical(kind=c_bool) :: sucess !< Flag indicating if the read was sucessful
end function open_and_parse_file_wrap

!> @brief Private c function that checks if a file_id is valid (see yaml_parser_binding.c)
!! @return Flag indicating if the file_id is valid
function is_valid_file_id(file_id) bind(c) &
   result(is_valid)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   logical(kind=c_bool) :: is_valid !< Flag indicating if the file_id is valid
end function is_valid_file_id

!> @brief Private c function that gets the number of key-value pairs in a block (see yaml_parser_binding.c)
!! @return Number of key-value pairs in this block
function get_nkeys_binding(file_id, block_id) bind(c) &
   result(nkeys)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer(kind=c_int), intent(in) :: block_id !< Id of the parent_block
   integer(kind=c_int) :: nkeys
end function get_nkeys_binding

!> @brief Private c function that gets the ids of the key-value pairs in a block (see yaml_parser_binding.c)
subroutine get_key_ids_binding(file_id, block_id, key_ids) bind(c)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer(kind=c_int), intent(in) :: block_id !< Id of the parent_block
   integer(kind=c_int), intent(inout) :: key_ids(*) !< Ids of the key-value pairs
end subroutine get_key_ids_binding

!> @brief Private c function that checks if a key_id is valid (see yaml_parser_binding.c)
!! @return Flag indicating if the key_id is valid
function is_valid_key_id(file_id, key_id) bind(c) &
   result(is_valid)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer(kind=c_int), intent(in) :: key_id !< Key id to check if valid
   logical(kind=c_bool) :: is_valid !< Flag indicating if the file_id is valid
end function is_valid_key_id

!> @brief Private c function that get the key from a key_id in a yaml file
!! @return Name of the key obtained
function get_key(file_id, key_id) bind(c) &
   result(key_name)
   use iso_c_binding, only: c_ptr, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer(kind=c_int), intent(in) :: key_id !< Id of the key-value pair of interest
   type(c_ptr) :: key_name
end function get_key

!> @brief Private c function that get the value from a key_id in a yaml file
!! @return String containing the value obtained
function get_value(file_id, key_id) bind(c) &
   result(key_value)
   use iso_c_binding, only: c_ptr, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer(kind=c_int), intent(in) :: key_id !< Id of the key-value pair of interest
   type(c_ptr) :: key_value
end function get_value

!> @brief Private c function that determines they value of a key in yaml_file (see yaml_parser_binding.c)
!! @return c pointer with the value obtained
function get_value_from_key_wrap(file_id, block_id, key_name, sucess) bind(c) &
   result(key_value2)

   use iso_c_binding, only: c_ptr, c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file to search
   integer(kind=c_int), intent(in) :: block_id !< ID corresponding to the block you want the key for
   character(kind=c_char), intent(in) :: key_name(*) !< Name of the key you want the value for
   integer(kind=c_int), intent(out) :: sucess !< Flag indicating if the call was sucessful
   type(c_ptr) :: key_value2
end function get_value_from_key_wrap

!> @brief Private c function that determines the number of blocks with block_name in the yaml file
!! (see yaml_parser_binding.c)
!! @return Number of blocks with block_name
function get_num_blocks_all(file_id, block_name) bind(c) &
   result(nblocks)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file to search
   character(kind=c_char), intent(in) :: block_name(*) !< The name of the block you are looking for

   integer(kind=c_int) :: nblocks
end function get_num_blocks_all

!> @brief Private c function that determines the number of blocks with block_name that belong to
!! a parent block with parent_block_id in the yaml file (see yaml_parser_binding.c)
!! @return Number of blocks with block_name
function get_num_blocks_child(file_id, block_name, parent_block_id) bind(c) &
   result(nblocks)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file to search
   character(kind=c_char), intent(in) :: block_name(*) !< The name of the block you are looking for
   integer(kind=c_int) :: parent_block_id !< Id of the parent block

   integer(kind=c_int) :: nblocks
end function get_num_blocks_child

!> @brief Private c function that gets the the ids of the blocks with block_name in the yaml file
!! (see yaml_parser_binding.c)
subroutine get_block_ids_all(file_id, block_name, block_ids) bind(c)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file to search
   character(kind=c_char), intent(in) :: block_name(*) !< The name of the block you are looking for
   integer(kind=c_int), intent(inout) :: block_ids(*) !< Id of the parent_block
end subroutine get_block_ids_all

!> @brief Private c function that gets the the ids of the blocks with block_name and that
!! belong to a parent block id in the yaml file (see yaml_parser_binding.c)
subroutine get_block_ids_child(file_id, block_name, block_ids, parent_block_id) bind(c)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file to search
   character(kind=c_char), intent(in) :: block_name(*) !< The name of the block you are looking for
   integer(kind=c_int), intent(inout) :: block_ids(*) !< Id of the parent_block
   integer(kind=c_int) :: parent_block_id !< Id of the parent block
end subroutine get_block_ids_child

!> @brief Private c function that checks if a block_id is valid (see yaml_parser_binding.c)
!! @return Flag indicating if the block_id is valid
function is_valid_block_id(file_id, block_id) bind(c) &
   result(is_valid)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer(kind=c_int), intent(in) :: block_id !< Block id to check if valid
   logical(kind=c_bool) :: is_valid !< Flag indicating if the file_id is valid
end function is_valid_block_id

end interface

!> @addtogroup yaml_parser_mod
!> @{
contains

!> @brief Opens and parses a yaml file
!! @return A file id corresponding to the file that was opened
function open_and_parse_file(filename) &
   result(file_id)

   character(len=*), intent(in) :: filename !< Filename of the yaml file
   logical :: sucess !< Flag indicating if the read was sucessful
   logical :: yaml_exists !< Flag indicating whether the yaml exists

   integer :: file_id

   inquire(file=trim(filename), EXIST=yaml_exists)
   if (.not. yaml_exists) then
      file_id = 999
      call mpp_error(NOTE, "The yaml file:"//trim(filename)//" does not exist, hopefully this is your intent!")
      return
   end if
   sucess = open_and_parse_file_wrap(trim(filename)//c_null_char, file_id)
   if (.not. sucess) call mpp_error(FATAL, "Error opening the yaml file:"//trim(filename)//". Check the file!")

end function open_and_parse_file

!> @brief Gets the key from a file id
subroutine get_key_name(file_id, key_id, key_name)
   integer, intent(in) :: key_id !< Id of the key-value pair of interest
   integer, intent(in) :: file_id !< File id of the yaml file to search
   character(len=*), intent(out) :: key_name

   if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
       &  "The file id in your get_key_name call is invalid! Check your call.")
   if (.not. is_valid_key_id(file_id, key_id)) call mpp_error(FATAL, &
       &  "The key id in your get_key_name call is invalid! Check your call.")

   key_name = fms_c2f_string(get_key(file_id, key_id))

end subroutine get_key_name

!> @brief Gets the value from a file id
subroutine get_key_value(file_id, key_id, key_value)
   integer, intent(in) :: key_id !< Id of the key-value pair of interest
   integer, intent(in) :: file_id !< File id of the yaml file to search
   character(len=*), intent(out) :: key_value

   if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
       &  "The file id in your get_key_value call is invalid! Check your call.")
   if (.not. is_valid_key_id(file_id, key_id)) call mpp_error(FATAL, &
       &  "The key id in your get_key_value call is invalid! Check your call.")

   key_value = fms_c2f_string(get_value(file_id, key_id))

end subroutine get_key_value

!> @brief Used to dermine the value of a key from a keyname
subroutine get_value_from_key_0d(file_id, block_id, key_name, key_value, is_optional)
   integer, intent(in) :: file_id !< File id of the yaml file to search
   integer, intent(in) :: block_id !< ID corresponding to the block you want the key for
   character(len=*), intent(in) :: key_name !< Name of the key you want the value for
   class(*), intent(inout):: key_value !< Value of the key
   logical, intent(in), optional :: is_optional !< Flag indicating if it is okay for they key to not exist.
                                                !! If the key does not exist key_value will not be set, so it
                                                !! is the user's responsibility to initialize it before the call

   character(len=255) :: buffer !< String buffer with the value

   type(c_ptr) :: c_buffer !< c pointer with the value
   integer(kind=c_int) :: sucess !< Flag indicating if the value was obtained sucessfully
   logical :: optional !< Flag indicating that the key was optional
   integer :: err_unit !< integer with io error

   optional = .false.
   if (present(is_optional)) optional = is_optional

   if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
       &  "The file id in your get_value_from_key call is invalid! Check your call.")
   if (.not. is_valid_block_id(file_id, block_id)) call mpp_error(FATAL, &
       &  "The block id in your get_value_from_key call is invalid! Check your call.")

   c_buffer = get_value_from_key_wrap(file_id, block_id, trim(key_name)//c_null_char, sucess)
   if (sucess == 1) then
     buffer = fms_c2f_string(c_buffer)

     select type (key_value)
       type is (integer(kind=i4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to i4")
       type is (integer(kind=i8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to i8")
       type is (real(kind=r4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to r4")
       type is (real(kind=r8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to r8")
       type is (character(len=*))
          call string_copy(key_value, buffer)
       type is (logical)
          if (lowercase(trim(buffer)) == "false") then
            key_value = .false.
          elseif (lowercase(trim(buffer)) == "true") then
            key_value = .true.
          else
            call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to logical")
          endif
     class default
       call mpp_error(FATAL, "The type of your buffer in your get_value_from_key call for key "//trim(key_name)//&
                            &" is not supported. Only i4, i8, r4, r8 and strings are supported.")
     end select
   else
     if(.not. optional) call mpp_error(FATAL, "Error getting the value for key:"//trim(key_name))
   endif

end subroutine get_value_from_key_0d

!> @brief Used' to dermine the 1D value of a key from a keyname
subroutine get_value_from_key_1d(file_id, block_id, key_name, key_value, is_optional)
   integer, intent(in) :: file_id !< File id of the yaml file to search
   integer, intent(in) :: block_id !< ID corresponding to the block you want the key for
   character(len=*), intent(in) :: key_name !< Name of the key you want the value for
   class(*), intent(inout):: key_value(:) !< Value of the key
   logical, intent(in), optional :: is_optional !< Flag indicating if it is okay for they key' to not exist.
                                                !! If the key does not exist key_value will not be set, so it
                                                !! is the user's responsibility to initialize it before the call

   character(len=255) :: buffer !< String buffer with the value

   type(c_ptr) :: c_buffer !< c pointer with the value
   integer(kind=c_int) :: sucess !< Flag indicating if the value was obtained sucessfully
   logical :: optional !< Flag indicating that the key was optional
   integer :: err_unit !< integer with io error

   optional=.false.
   if (present(is_optional)) optional = is_optional

   if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
       &  "The file id in your get_value_from_key call is invalid! Check your call.")
   if (.not. is_valid_block_id(file_id, block_id)) call mpp_error(FATAL, &
       &  "The block id in your get_value_from_key call is invalid! Check your call.")

   c_buffer = get_value_from_key_wrap(file_id, block_id, trim(key_name)//c_null_char, sucess)
   if (sucess == 1) then
     buffer = fms_c2f_string(c_buffer)

     select type (key_value)
       type is (integer(kind=i4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to i4")
       type is (integer(kind=i8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to i8")
       type is (real(kind=r4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to r4")
       type is (real(kind=r8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"// &
              & trim(key_name)//" Error converting '"//trim(buffer)//"' to r8")
       type is (character(len=*))
          call mpp_error(FATAL, "get_value_from_key 1d string variables are not supported. Contact developers")
     class default
       call mpp_error(FATAL, "The type of your buffer in your get_value_from_key call for key "//trim(key_name)//&
                            &" is not supported. Only i4, i8, r4, r8 and strings are supported.")
     end select
   else
     if(.not. optional) call mpp_error(FATAL, "Error getting the value for key:"//trim(key_name))
   endif
end subroutine get_value_from_key_1d

!> @brief Determines the number of blocks with block_name in the yaml file
!! If parent_block_id is present, it only counts those that belong to that block
!! @return Number of blocks with block_name
function get_num_blocks(file_id, block_name, parent_block_id) &
    result(nblocks)

    integer, intent(in) :: file_id !< File id of the yaml file to search
    character(len=*), intent(in) :: block_name !< The name of the block you are looking for
    integer, intent(in), optional :: parent_block_id !< Id of the parent block
    integer :: nblocks

    if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
        &  "The file id in your get_num_blocks call is invalid! Check your call.")

    if (.not. present(parent_block_id)) then
       nblocks=get_num_blocks_all(file_id, trim(block_name)//c_null_char)
    else
       if (.not. is_valid_block_id(file_id, parent_block_id)) call mpp_error(FATAL, &
           &  "The parent_block id in your get_num_blocks call is invalid! Check your call.")
       nblocks=get_num_blocks_child(file_id, trim(block_name)//c_null_char, parent_block_id)
    endif
end function get_num_blocks

!> @brief Gets the the ids of the blocks with block_name in the yaml file
!! If parent_block_id is present, it only gets those that belong to that block
subroutine get_block_ids(file_id, block_name, block_ids, parent_block_id)

    integer, intent(in) :: file_id !< File id of the yaml file to search
    character(len=*), intent(in) :: block_name !< The name of the block you are looking for
    integer, intent(inout) :: block_ids(:) !< Id of blocks with block_name
    integer, intent(in), optional :: parent_block_id !< Id of the parent_block
    integer :: nblocks_id
    integer :: nblocks

    if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
        &  "The file id in your get_block_ids call is invalid! Check your call.")

    nblocks_id = size(block_ids)
    nblocks = get_num_blocks(file_id, block_name, parent_block_id)
    if (nblocks .ne. nblocks_id) call mpp_error(FATAL, "The size of your block_ids array is not correct")

    if (.not. present(parent_block_id)) then
       call get_block_ids_all(file_id, trim(block_name)//c_null_char, block_ids)
    else
       if (.not. is_valid_block_id(file_id, parent_block_id)) call mpp_error(FATAL, &
           &  "The parent_block id in your get_block_ids call is invalid! Check your call.")
       call get_block_ids_child(file_id, trim(block_name)//c_null_char, block_ids, parent_block_id)
    endif
end subroutine get_block_ids

!> @brief Gets the number of key-value pairs in a block
!! @return Number of key-value pairs in this block
function get_nkeys(file_id, block_id) &
   result(nkeys)
   integer, intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer, intent(in) :: block_id !< Id of the parent_block
   integer :: nkeys

    if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
        &  "The file id in your get_nkeys call is invalid! Check your call.")
    if (.not. is_valid_block_id(file_id, block_id)) call mpp_error(FATAL, &
        &  "The block id in your get_nkeys call is invalid! Check your call.")

    nkeys = get_nkeys_binding(file_id, block_id)
end function get_nkeys

!> @brief Gets the ids of the key-value pairs in a block
subroutine get_key_ids (file_id, block_id, key_ids)
   integer, intent(in) :: file_id !< File id corresponding to the yaml file that was opened
   integer, intent(in) :: block_id !< Id of the parent_block
   integer, intent(inout) :: key_ids(:) !< Ids of the key-value pairs

   integer :: nkey_ids !< Size of key_ids
   integer :: nkeys !< Actual number of keys

   if (.not. is_valid_file_id(file_id)) call mpp_error(FATAL, &
       &  "The file id in your get_key_ids call is invalid! Check your call.")
   if (.not. is_valid_block_id(file_id, block_id)) call mpp_error(FATAL, &
       &  "The block id in your get_key_ids call is invalid! Check your call.")

   nkey_ids = size(key_ids)
   nkeys = get_nkeys(file_id, block_id)

   if (nkeys .ne. nkey_ids) call mpp_error(FATAL, "The size of your key_ids array is not correct.")

   call get_key_ids_binding (file_id, block_id, key_ids)
end subroutine get_key_ids

#endif
end module yaml_parser_mod
!> @}
! close documentation grouping
