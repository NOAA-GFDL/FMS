## Check is compiler supports "type is (character(len=*))"
AC_DEFUN([TYPE_IS_CHECK],[
AC_LANG_PUSH([Fortran])
AC_COMPILE_IFELSE([subroutine test_sub(ctype)
  class(*), intent(out) :: ctype

  select type(ctype)
    type is (character(len=*))
      ctype(:) = ""
  end select
end subroutine test_sub], [], AC_MSG_ERROR([This compiler doesn't support "type is (character(len=*)")]))

AC_LANG_POP([Fortran])
])
