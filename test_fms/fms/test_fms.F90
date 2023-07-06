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
 use mpp_mod, only : mpp_error, fatal, note, mpp_init, stderr
 use fms_mod, only : fms_init, fms_end
 use fms_mod, only : fms_c2f_string
 use fms_mod, only : fms_cstring2cpointer
 use fms_mod, only : monotonic_array
 use platform_mod, only : r4_kind, r8_kind
 use fms_string_utils_mod, only : string, stringify
 use test_fms_mod
 use, intrinsic :: iso_c_binding

 integer :: i !< Integer
 character(len=16) :: answer !< expected answer
 character(len=16) :: test !< Test string
 character(kind=c_char) :: Cstring (17)!< C string to convert
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
!!!!!!!!!!!!!!!!!!!!
! Test the c string to c pointer conversion
 test = "                "
 answer = '100'
 Cstring =  "                "
 Cstring(1) = "1"
 Cstring(2) = "0"
 Cstring(3) = "0"
 Cstring(17) = c_null_char
 call mpp_error(NOTE,"Testing fms_cstring2cpointer and fms_c2f_string")
! test = fms_c2f_string(fms_cstring2cpointer(c_char_"100             "//c_null_char))
 test = fms_c2f_string(Cstring)
 if (trim(answer) .eq. trim(test)) then
         call mpp_error(NOTE, trim(test)//" matches "//trim(answer))
 else
         call mpp_error(FATAL, trim(test)//" does not match "//trim(answer))
 endif

 call test_monotonic_array_r4
 call test_monotonic_array_r8

 call fms_end()

contains


#include "test_fms_r4.fh"
#include "test_fms_r8.fh"

end program test_fms
