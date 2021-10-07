module test_fms_mod
 use, intrinsic :: iso_c_binding

        interface 
                function strPoint () bind(C,name="strPoint") result(Cstring)
                   use, intrinsic :: iso_c_binding
                        type (c_ptr) :: Cstring !< String returned from this function "100"
                end function strPoint
        end interface

end module test_fms_mod

program test_fms
 use mpp_mod, only : mpp_error, fatal, note, mpp_init
 use fms_mod, only : fms_init, string, fms_end
 use fms_mod, only : fms_c2f_string
 use test_fms_mod
 use, intrinsic :: iso_c_binding

 integer :: i !< Integer
 character(len=16) :: answer !< expected answer
 character(len=16) :: test !< Test string
 character(len=:,kind=c_char), pointer :: Cstring !< C string to convert
 type(c_ptr), pointer :: Cptr !< C pointer to string

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

!Test the C string to F string function
 test = "                "
 answer = '100'
! Cptr = strPoint()
 test = fms_c2f_string(strPoint())

 write (6,*) "Cstring test = ", trim(test)
 if (trim(answer) .eq. trim(test)) then
         call mpp_error(NOTE, trim(test)//" matches "//trim(answer))
 else
         call mpp_error(FATAL, trim(test)//" does not match "//trim(answer))
 endif

 call fms_end()
 
end program test_fms
