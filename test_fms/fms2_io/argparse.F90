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

!> @brief A simple command line argument parsing module for Fortran programs.
module argparse
use, intrinsic :: iso_fortran_env, only: error_unit
implicit none
private


!> @brief Argument type.
type, private :: Argument_t
  character(len=32) :: short_name !< Short name (must start with a
                                  !! - if it is optional).
  character(len=32) :: long_name !< Long name (must start witha
                                 !! --).
  character(len=512) :: description !< Description that will appear
                                    !! in the help message.
  character(len=512) :: val !< Buffer where the input value is stored.
  logical :: positional !< Flag telling if the argument is positional.
  logical :: requires_val !< Flag telling if the argument requires a
                          !! value.
  type(Argument_t), pointer :: head !< Makes this type a linked list.
endtype Argument_t


!> @brief Parser type.
type, public :: Parser_t
  type(Argument_t),pointer :: args !< Linked list of arguments.
                                   !! Will be made a dictionary
                                   !! in the future.
  integer :: argc !< Total number of arguments.
  integer :: num_pos_args !< Number of positional arguments.
endtype Parser_t


public :: get_parser
public :: destroy_parser
public :: add_argument
public :: parse_args
public :: get_argument


contains


!> @brief Initialize a parser type.
!! @return An initialized parser type with a help option.
function get_parser() result(p)

  type(Parser_t) :: p

  p%args => null()
  p%argc = 0
  p%num_pos_args = 0
  call add_argument(p, "-h", "Print this help message", .false., "--help")
end function get_parser

!> @brief Release memory used by a parser type.
subroutine destroy_parser(parser)

  type(Parser_t), intent(inout) :: parser !< Parser object.

  type(Argument_t), pointer :: p
  type(Argument_t), pointer :: p2

  p => parser%args
  do while (associated(p))
    p2 => p%head
    deallocate(p)
    p => p2
  enddo
  parser%args => null()
  parser%argc = 0
  parser%num_pos_args = 0
end subroutine destroy_parser

!> @brief Add an argument to the parser object.
subroutine add_argument(parser,short_name,description,requires_val, long_name)

  type(Parser_t), intent(inout) :: parser !< Parser object.
  character(len=*), intent(in) :: short_name !< Short name.
  character(len=*), intent(in) :: description !< What the option does.
  logical, intent(in), optional :: requires_val !< Does the option require
                                                !! a value?
  character(len=*), intent(in), optional :: long_name !< Long name.
  type(Argument_t), pointer :: arg
  type(Argument_t), pointer :: a

  !Create a new arg linked list node.
  allocate(arg)
  arg%head => null()
  call string_copy(arg%short_name, short_name)
  call string_copy(arg%description, description)
  parser%argc = parser%argc + 1
  arg%long_name = ""
  arg%val = ""

  !Optional arguments must have short names that begin with a -, or
  !else they are assumed to be positional arguments.
  if (short_name(1:1) .eq. "-") then
    arg%positional = .false.
    if (present(long_name)) then

      !Long names are required to start with --.
      if (long_name(1:2) .ne. "--") then
        call error("argument long names must start with --")
      endif
      call string_copy(arg%long_name, long_name)
    endif
    if (.not. present(requires_val)) then
      call error("you must specify whether or not a value is" &
                 //" required for optional arguments.")
    endif
    arg%requires_val = requires_val
  else
    arg%positional = .true.
    arg%requires_val = .false.
    parser%num_pos_args = parser%num_pos_args + 1
  endif

  !Add the new argument to the end of the linked list.  Trap any
  !duplicate arguments (this results in an error).
  !To do: switch from a linked list to a dictionary to increase
  !efficiency.
  if (associated(parser%args)) then
    a => parser%args
    do while (.true.)
      call catch_duplicate_args(a, arg)
      if (associated(a%head)) then
        a => a%head
      else
        exit
      endif
    enddo
    a%head => arg
  else
    parser%args => arg
  endif
end subroutine add_argument


!> @brief Print a help message.
subroutine help_mesg(parser)

  type(Parser_t), intent(inout) :: parser !< Parser object.

  type(Argument_t), pointer :: arg

  !Write out positional arguments first.
  if (parser%num_pos_args .gt. 0) then
    write(error_unit,*) "Positional arguments:"
    arg => parser%args
    do while (associated(arg))
      if (arg%positional) then
        call print_argument(arg)
      endif
      arg => arg%head
    enddo
    write(error_unit,*)
  endif

  !Write out optional arguments.
  if (parser%argc .gt. parser%num_pos_args) then
    write(error_unit,*) "Optional arguments:"
    arg => parser%args
    do while (associated(arg))
      if (.not. arg%positional) then
        call print_argument(arg)
      endif
      arg => arg%head
    enddo
    write(error_unit,*)
  endif
end subroutine help_mesg


!> @brief Print a usage message.
subroutine usage_mesg(parser)

  type(Parser_t),intent(inout) :: parser !< Parser object.

  type(Argument_t), pointer :: arg
  character(len=512) :: buf
  character(len=512) :: b

  call get_command_argument(0, buf)
  call string_copy(buf,"Usage: "//trim(buf))
  if (parser%argc .gt. 0) then
    arg => parser%args
    do while (associated(arg))
      if (arg%positional) then
        call string_copy(b, arg%short_name)
      else
        !Print optional arguments as [-shortname|--longname <value>].
        call string_copy(b, "["//trim(arg%short_name))
        if (len_trim(arg%long_name) .gt. 0) then
          call string_copy(b, trim(b)//"|"//trim(arg%long_name))
        endif
        if (arg%requires_val) then
          call string_copy(b, trim(b)//" <value>")
        endif
        call string_copy(b, trim(b)//"]")
      endif
      call string_copy(buf, trim(buf)//" "//trim(b))
      arg => arg%head
    enddo
  endif
  write(error_unit, '(a)') trim(buf)
end subroutine usage_mesg


!> @brief Loop through the supplied command line arguments and store
!! them in the parser object's arguments linked list.
subroutine parse_args(parser)

  type(Parser_t), intent(inout) :: parser !< Parser object.

  character(len=128), dimension(:), allocatable :: args
  type(Argument_t), pointer :: arg
  integer :: argc
  integer :: i
  logical :: arg_found
  logical :: skip
  integer :: cl_pos_arg_counter
  integer :: cl_opt_arg_counter

  cl_pos_arg_counter = 0
  cl_opt_arg_counter = parser%num_pos_args

  argc = command_argument_count()
  allocate(args(argc))
  do i = 1, argc
    call get_command_argument(i, args(i))
  enddo

  arg_found = .false.
  skip = .false.
  do i = 1, argc
    if (skip) then
      arg_found = .false.
      skip = .false.
      cycle
    endif
    if (trim(args(i)) .eq. "-h" .or. trim(args(i)) .eq. "--help") then
      call usage_mesg(parser)
      call help_mesg(parser)
      stop
    endif
    arg => parser%args
    do while (associated(arg))
      if (.not. arg%positional) then
        if (trim(arg%long_name) .eq. trim(args(i))) then
          call string_copy(args(i), arg%short_name)
        endif
        if (trim(arg%short_name) .eq. trim(args(i))) then
          cl_opt_arg_counter = cl_opt_arg_counter + 1
          if (cl_opt_arg_counter .gt. parser%argc) then
            call usage_mesg(parser)
            call error("too many input arguments.")
          endif
          if (trim(arg%val) .ne. "") then
            call error("duplicate arguments given.")
          endif
          if (arg%requires_val) then
            if (i .eq. argc) then
              call usage_mesg(parser)
              call error("too few input arguments.")
            endif
            skip = .true.
            call string_copy(arg%val, args(i+1))
          else
            call string_copy(arg%val, "present")
          endif
          arg_found = .true.
          exit
        endif
      endif
      arg => arg%head
    enddo
    if (arg_found) then
      arg_found = .false.
      cycle
    endif
    cl_pos_arg_counter = cl_pos_arg_counter + 1
    if (cl_pos_arg_counter .gt. parser%num_pos_args) then
      call error("too many positional arguments.")
    endif
    arg => parser%args
    do while (associated(arg))
      if (arg%positional .and. trim(arg%val) .eq. "") then
        call string_copy(arg%val,args(i))
        exit
      endif
      arg => arg%head
    enddo
  enddo
  if (cl_pos_arg_counter .ne. parser%num_pos_args) then
    call usage_mesg(parser)
    call error("incorrect number of positional arguments.")
  endif
  deallocate(args)
end subroutine parse_args


!> @brief Given an argument's shortname, return the value that was passed on
!!        the command line.  If the argument is optional and not,
!!        not supplied, return the string "not present".  Note that this
!!        routine must be called after parse_args.
subroutine get_argument(parser, short_name, val)

  type(Parser_t), intent(in) :: parser
  character(len=*), intent(in) :: short_name
  character(len=*), intent(inout) :: val

  type(Argument_t), pointer :: arg

  arg => parser%args
  do while (associated(arg))
    if (trim(arg%short_name) .eq. trim(short_name)) then
      if (trim(arg%val) .eq. "") then
        call string_copy(val, "not present")
      else
        call string_copy(val, arg%val)
      endif
      return
    endif
    arg => arg%head
  enddo
  call error("unrecognized argument "//trim(short_name)//".")
end subroutine get_argument


!> @brief Print an argument name and its description to stderr.
subroutine print_argument(arg)

  type(Argument_t), intent(in) :: arg !< Argument.

  character(len=32) :: buf

  call string_copy(buf, arg%short_name)
  if (len_trim(arg%long_name) .gt. 0) then
    call string_copy(buf, trim(buf)//","//arg%long_name)
  endif
  if (arg%requires_val) then
    write(error_unit, '(a32,1x,a5,7x,a)') trim(buf),"value", trim(arg%description)
  else
    write(error_unit, '(a32,13x,a)') trim(buf), trim(arg%description)
  endif
end subroutine print_argument


!> @brief Throw an error if the two input arguments have the same
!! short and/or long names.
subroutine catch_duplicate_args(arg1, arg2)

  type(Argument_t), intent(in) :: arg1 !< First argument.
  type(Argument_t), intent(in) :: arg2 !< Second argument.

  if (trim(arg1%short_name) .eq. trim(arg2%short_name)) then
    call error("argument "//trim(arg1%short_name)//" already exists.")
  endif
  if (len_trim(arg1%long_name) .ge. 1 .and. trim(arg1%long_name) .eq. &
      trim(arg2%long_name)) then
    call error("argument "//trim(arg1%long_name)//" already exists.")
  endif
end subroutine catch_duplicate_args


!> @brief Print a message to stderr, then stop the program.
subroutine error(mesg)

  character(len=*),intent(in) :: mesg !< Message that will be printed to
                                      !! stderr.

  write(error_unit,*) "Error: "//trim(mesg)
  stop 1
end subroutine error


!> @brief Safely copy a string from one buffer to another.
subroutine string_copy(dest, source)

  character(len=*), intent(inout) :: dest !< Destination string.
  character(len=*), intent(in) :: source !< Source string.

  if (len_trim(source) .gt. len(dest)) then
    call error("The input destination string is not big enough to" &
               //" to hold the input source string.")
  endif
  dest = ""
  dest = trim(source)
end subroutine string_copy


end module argparse


#ifdef UNITTEST
program test
use, intrinsic :: iso_fortran_env, only: output_unit
use argparse
implicit none

type(Parser_t) :: parser

parser = get_parser()
call add_argument(parser, "pos1", "first positional arg.")
call add_argument(parser, "pos2", "second positional arg.")
call add_argument(parser, "-v", "Verbose mode", requires_val=.false., &
                  long_name="--verbose")
call add_argument(parser, "-a", "Optional arg", requires_val=.true.)
call parse_args(parser)
call show_arg(parser, "pos1")
call show_arg(parser, "pos2")
call show_arg(parser, "-v")
call show_arg(parser, "-a")


contains


subroutine show_arg(parser, shortname)

  type(Parser_t), intent(in) :: parser
  character(len=*), intent(in) :: shortname

  character(len=128) :: buf

  call get_argument(parser, shortname, buf)
  write(output_unit, *) trim(shortname)//" = "//trim(buf)
end subroutine show_arg


end program test
#endif
