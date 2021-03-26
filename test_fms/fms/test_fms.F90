program test_fms
 use mpp_mod, only : mpp_error, fatal, note, mpp_init
 use fms_mod, only : fms_init, string

 integer :: i !< Integer
 character(len=16) :: answer !< expected answer
 character(len=16) :: test !< Test string

 call mpp_init()
 call fms_init()

 test = "                "
 answer = "          "
 answer = '100'
 i = 100
 test = string(i)
 write (6,*) "Integer test = ", trim(test)
 if (trim(answer) .eq. trim(test)) then
         call mpp_error(NOTE, trim(test)//" matches "//trim(answer))
 else
         call mpp_error(FATAL, trim(test)//" does not match "//trim(answer))
 endif
 
end program test_fms
