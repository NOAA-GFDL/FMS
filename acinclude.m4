dnl Some autoconf macros to help build GFDL Fortran projects.

dnl Ed Hartnett, 9/11/19

dnl Check that Fortran compiler can hanlde long lines of code.
AC_DEFUN([GFDL_CHECK_FORTRAN_LONG_LINES],
[
# Tell autoconf that language tests are in fortran.
AC_LANG_PUSH(Fortran)

AC_FC_SRCEXT(F90)
AC_MSG_CHECKING([whether Fortran compiler can handle long lines of code])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
        subroutine foo(bar)
        integer, intent(in) :: bar
        print *,'7890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890'
        end subroutine foo
        ])],
        [test_result=yes],
        [test_result=no])
AC_MSG_RESULT($test_result)
if test $test_result = no; then
   AC_MSG_ERROR([Fortran compiler must be able to handle long lines of code. \
   Set FCFLAGS to allow this (-ffree-line-length-none for gfortran)])
fi

# Done testing Fortran compiler.
AC_LANG_POP(Fortran)
])
