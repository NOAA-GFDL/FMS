module yaml_parser_mod

implicit none

public :: do_stuff
public :: read_and_parse_file

interface
function read_and_parse_file(filename) bind(c) &
   result(sucess)
   use iso_c_binding, only: c_char, c_int, c_bool
   character(kind=c_char), intent(in) :: filename(*)
   logical(kind=c_bool) :: sucess
end function read_and_parse_file
end interface

contains

subroutine do_stuff()
   print *, "Doing stuff"

#ifdef use_yaml
   print *, "Very important stuff"
#endif

end subroutine do_stuff

end module yaml_parser_mod
