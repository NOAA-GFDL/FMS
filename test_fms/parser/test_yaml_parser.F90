program test_yaml_parser

use yaml_parser_mod
use mpp_mod
use fms_mod, only : fms_init, fms_end

implicit none

call do_stuff()

#ifdef use_yaml

call fms_init()
if (read_and_parse_file("diag_table.yaml")) then
    print *, "The yaml file was read bro ^"
else
   call mpp_error(FATAL, "The file was not opened sucessfully, get help!")
endif
call fms_end()

#endif

end program test_yaml_parser
