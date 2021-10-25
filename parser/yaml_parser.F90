module yaml_parser_mod

#ifdef use_yaml
use fms_mod, only: fms_c2f_string
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

!> @brief Dermine the value of a key from a keyname
interface get_value_from_key
  module procedure get_value_from_key_0d
  module procedure get_value_from_key_1d
end interface get_value_from_key

interface

!> @brief Private c function that opens and parses a yaml file (see yaml_parser_binding.c)
!! @return Flag indicating if the read was sucessful
function open_and_parse_file_wrap(filename, file_id) bind(c) &
   result(sucess)
   use iso_c_binding, only: c_char, c_int, c_bool
   character(kind=c_char), intent(in) :: filename(*) !< Filename of the yaml file
   integer(kind=c_int), intent(out) :: file_id !< File id corresponding' to the yaml file that was opened
   logical(kind=c_bool) :: sucess !< Flag indicating if the read was sucessful
end function open_and_parse_file_wrap

!> @brief c function that gets the number of key-value pairs in a block (see yaml_parser_binding.c)
!! @return Number of key-value pairs in this block
function get_nkeys(file_id, block_id) bind(c) &
   result(nkeys)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding' to the yaml file that was opened
   integer(kind=c_int), intent(in) :: block_id !< Id of the parent_block
   integer(kind=c_int) :: nkeys
end function get_nkeys

!> @brief c function that gets the ids of the key-value pairs in a block (see yaml_parser_binding.c)
subroutine get_key_ids(file_id, block_id, key_ids) bind(c)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding' to the yaml file that was opened
   integer(kind=c_int), intent(in) :: block_id !< Id of the parent_block
   integer(kind=c_int), intent(inout) :: key_ids(*) !< Ids of the key-value pairs
end subroutine get_key_ids

!> @brief Private c function that get the key from a key_id in a yaml file
!! @return Name of the key obtained
function get_key(file_id, key_id) bind(c) &
   result(key_name)
   use iso_c_binding, only: c_ptr, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding' to the yaml file that was opened
   integer(kind=c_int), intent(in) :: key_id !< Id of the key-value pair of interest
   type(c_ptr) :: key_name
end function get_key

!> @brief Private c function that get the value from a key_id in a yaml file
!! @return String containing the value obtained
function get_value(file_id, key_id) bind(c) &
   result(key_value)
   use iso_c_binding, only: c_ptr, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id corresponding' to the yaml file that was opened
   integer(kind=c_int), intent(in) :: key_id !< Id of the key-value pair of interest
   type(c_ptr) :: key_value
end function get_value

!> @brief Private c function that determines they value of a key in yaml_file (see yaml_parser_binding.c)
!! @return c pointer with the value obtained
function get_value_from_key_wrap(file_id, block_id, key_name, sucess) bind(c) &
   result(key_value2)

   use iso_c_binding, only: c_ptr, c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file' to search
   integer(kind=c_int), intent(in) :: block_id !< ID corresponding' to the block you want the key for
   character(kind=c_char), intent(in) :: key_name !< Name of the key you want the value for
   logical(kind=c_bool), intent(out) :: sucess !< Flag indicating if the call was sucessful
   type(c_ptr) :: key_value2
end function get_value_from_key_wrap

!> @brief Private c function that determines the number of blocks with block_name in the yaml file
!! (see yaml_parser_binding.c)
!! @return Number of blocks with block_name
function get_num_blocks_all(file_id, block_name) bind(c) &
   result(nblocks)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file' to search
   character(kind=c_char), intent(in) :: block_name !< The name of the block you are looking for

   integer(kind=c_int) :: nblocks
end function get_num_blocks_all

!> @brief Private c function that determines the number of blocks with block_name that belong to
!! a parent block with parent_block_id in the yaml file (see yaml_parser_binding.c)
!! @return Number of blocks with block_name
function get_num_blocks_child(file_id, block_name, parent_block_id) bind(c) &
   result(nblocks)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file' to search
   character(kind=c_char), intent(in) :: block_name !< The name of the block you are looking for
   integer(kind=c_int) :: parent_block_id !< Id of the parent block

   integer(kind=c_int) :: nblocks
end function get_num_blocks_child

!> @brief Private c function that gets the the ids of the blocks with block_name in the yaml file
!! (see yaml_parser_binding.c)
subroutine get_block_ids_all(file_id, block_name, block_ids) bind(c)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file' to search
   character(kind=c_char), intent(in) :: block_name !< The name of the block you are looking for
   integer(kind=c_int), intent(inout) :: block_ids(*) !< Id of the parent_block
end subroutine get_block_ids_all

!> @brief Private c function that gets the the ids of the blocks with block_name and that
!! belong' to a parent block id in the yaml file (see yaml_parser_binding.c)
subroutine get_block_ids_child(file_id, block_name, block_ids, parent_block_id) bind(c)
   use iso_c_binding, only: c_char, c_int, c_bool
   integer(kind=c_int), intent(in) :: file_id !< File id of the yaml file' to search
   character(kind=c_char), intent(in) :: block_name !< The name of the block you are looking for
   integer(kind=c_int), intent(inout) :: block_ids(*) !< Id of the parent_block
   integer(kind=c_int) :: parent_block_id !< Id of the parent block
end subroutine get_block_ids_child

end interface

contains

!> @brief Opens and parses a yaml file
!! @return A file id corresponding' to the file that was opened
function open_and_parse_file(filename) &
   result(file_id)

   character(len=*), intent(in) :: filename !< Filename of the yaml file
   logical :: sucess !< Flag indicating if the read was sucessful

   integer :: file_id

   sucess = open_and_parse_file_wrap(filename, file_id)
   if (.not. sucess) call mpp_error(FATAL, "Error opening the yaml file:"//trim(filename)//". Check the file!")

end function open_and_parse_file

!> @brief Gets the key from a file id
subroutine get_key_name(file_id, key_id, key_name)
   integer, intent(in) :: key_id !< Id of the key-value pair of interest
   integer, intent(in) :: file_id !< File id of the yaml file' to search
   character(len=*), intent(out) :: key_name

   key_name = fms_c2f_string(get_key(file_id, key_id))

end subroutine get_key_name

!> @brief Gets the value from a file id
subroutine get_key_value(file_id, key_id, key_value)
   integer, intent(in) :: key_id !< Id of the key-value pair of interest
   integer, intent(in) :: file_id !< File id of the yaml file' to search
   character(len=*), intent(out) :: key_value

   key_value = fms_c2f_string(get_value(file_id, key_id))

end subroutine get_key_value

!> @brief Used to dermine the value of a key from a keyname
subroutine get_value_from_key_0d(file_id, block_id, key_name, key_value, is_optional)
   integer, intent(in) :: file_id !< File id of the yaml file' to search
   integer, intent(in) :: block_id !< ID corresponding' to the block you want the key for
   character(len=*), intent(in) :: key_name !< Name of the key you want the value for
   class(*), intent(inout):: key_value !< Value of the key
   logical, intent(in), optional :: is_optional !< Flag indicating if it is okay for they key' to not exist.
                                                !! If the key does not exist key_value will not be set, so it 
                                                !! is the user's responsibility' to initialize it before the call

   character(len=255) :: buffer !< String buffer with the value

   type(c_ptr) :: c_buffer !< c pointer with the value
   logical(kind=c_bool) :: sucess !< Flag indicating if the value was obtained sucessfully
   logical :: optional !< Flag indicating that the key was optional
   integer :: err_unit !< integer with io error

   if (present(is_optional)) optional = is_optional

   c_buffer = get_value_from_key_wrap(file_id, block_id, key_name, sucess)
   if (sucess) then
     buffer = fms_c2f_string(c_buffer)

     select type (key_value)
       type is (integer(kind=i4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to i4")
       type is (integer(kind=i8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to i8")
       type is (real(kind=r4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to r4")
       type is (real(kind=r8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to r8")
       type is (character(len=*))
          key_value = buffer
     class default
       call mpp_error(FATAL, "The type of your buffer in your get_value_from_key call for key "//trim(key_name)//&
                            &" is not supported. Only i4, i8, r4, r8 and strings are supported.")
     end select
   endif
   if(.not. sucess .and. .not. optional) call mpp_error(FATAL, "Error getting the value for key:"//trim(key_name))

end subroutine get_value_from_key_0d

!> @brief Used' to dermine the 1D value of a key from a keyname
subroutine get_value_from_key_1d(file_id, block_id, key_name, key_value, is_optional)
   integer, intent(in) :: file_id !< File id of the yaml file' to search
   integer, intent(in) :: block_id !< ID corresponding' to the block you want the key for
   character(len=*), intent(in) :: key_name !< Name of the key you want the value for
   class(*), intent(inout):: key_value(:) !< Value of the key
   logical, intent(in), optional :: is_optional !< Flag indicating if it is okay for they key' to not exist.
                                                !! If the key does not exist key_value will not be set, so it 
                                                !! is the user's responsibility' to initialize it before the call

   character(len=255) :: buffer !< String buffer with the value

   type(c_ptr) :: c_buffer !< c pointer with the value
   logical(kind=c_bool) :: sucess !< Flag indicating if the value was obtained sucessfully
   logical :: optional !< Flag indicating that the key was optional
   integer :: err_unit !< integer with io error

   if (present(is_optional)) optional = is_optional

   c_buffer = get_value_from_key_wrap(file_id, block_id, key_name, sucess)
   if (sucess) then
     buffer = fms_c2f_string(c_buffer)

     select type (key_value)
       type is (integer(kind=i4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to i4")
       type is (integer(kind=i8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to i8")
       type is (real(kind=r4_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to r4")
       type is (real(kind=r8_kind))
          read(buffer,*, iostat=err_unit) key_value
          if (err_unit .ne. 0) call mpp_error(FATAL, "Key:"//trim(key_name)//" Error converting '"//trim(buffer)//"' to r8")
       type is (character(len=*))
          call mpp_error(FATAL, "get_value_from_key 1d string variables are not supported. Contact developers")
     class default
       call mpp_error(FATAL, "The type of your buffer in your get_value_from_key call for key "//trim(key_name)//&
                            &" is not supported. Only i4, i8, r4, r8 and strings are supported.")
     end select
   endif
   if(.not. sucess .and. .not. optional) call mpp_error(FATAL, "Error getting the value for key:"//trim(key_name))

end subroutine get_value_from_key_1d

!> @brief Determines the number of blocks with block_name in the yaml file
!! If parent_block_id is present, it only counts those that belong' to that block
!! @return Number of blocks with block_name
function get_num_blocks(file_id, block_name, parent_block_id) &
    result(nblocks)

    integer, intent(in) :: file_id !< File id of the yaml file' to search
    character(len=*), intent(in) :: block_name !< The name of the block you are looking for
    integer, intent(in), optional :: parent_block_id !< Id of the parent block
    integer :: nblocks

    if (.not. present(parent_block_id)) then
       nblocks=get_num_blocks_all(file_id, block_name)
    else
       nblocks=get_num_blocks_child(file_id, block_name, parent_block_id)
    endif
end function get_num_blocks

!> @brief Gets the the ids of the blocks with block_name in the yaml file
!! If parent_block_id is present, it only gets those that belong' to that block
subroutine get_block_ids(file_id, block_name, block_ids, parent_block_id)

    integer, intent(in) :: file_id !< File id of the yaml file' to search
    character(len=*), intent(in) :: block_name !< The name of the block you are looking for
    integer, intent(inout) :: block_ids(:) !< Id of blocks with block_name
    integer, intent(in), optional :: parent_block_id !< Id of the parent_block

    if (.not. present(parent_block_id)) then
       call get_block_ids_all(file_id, block_name, block_ids)
    else
       call get_block_ids_child(file_id, block_name, block_ids, parent_block_id)
    endif
end subroutine get_block_ids

#endif
end module yaml_parser_mod
