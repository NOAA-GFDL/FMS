dnl Some autoconf macros to help build GFDL Fortran projects.


dnl Check that Fortran compiler can hanlde long lines of code.
dnl Ed Hartnett, 9/11/19
AC_DEFUN([GFDL_CHECK_FORTRAN_LONG_LINES],
[
# Tell autoconf that language tests are in fortran.
AC_LANG_PUSH(Fortran)

# Check a fortran compile with a test program with this file
# extension.
AC_FC_SRCEXT(F90)

# Attempt to compile test program with 160 char line of code.
AC_MSG_CHECKING([whether Fortran compiler can handle long lines of code])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
        subroutine foo(bar)
        integer, intent(in) :: bar
        print *,'789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890'
        end subroutine foo
        ])],
        [test_result=yes],
        [test_result=no])
AC_MSG_RESULT($test_result)

# If it failed, error out of configure with a helpful message.
if test $test_result = no; then
   AC_MSG_ERROR([Fortran compiler must be able to handle long lines of code. \
   Set FCFLAGS to allow this (-ffree-line-length-none for gfortran)])
fi

# Done testing Fortran compiler.
AC_LANG_POP(Fortran)
])


dnl Check that Fortran compiler is using 8-byte reals.
dnl Ed Hartnett, 9/11/19
AC_DEFUN([GFDL_CHECK_FORTRAN_REAL_8],
[
# Tell autoconf that language tests are in fortran.
AC_LANG_PUSH(Fortran)

# Check a fortran compile with a test program with this file
# extension.
AC_FC_SRCEXT(F90)

# Check the size of real.
AC_MSG_CHECKING([whether Fortran compiler has 8-byte real])
AC_RUN_IFELSE([AC_LANG_SOURCE([
        program foo
        real var
        if (sizeof(var) .ne. 8) then
        stop 2
        end if
        end program foo
        ])],
        [test_result=yes],
        [test_result=no],
        [AC_MSG_WARN([Test for 8-byte real cannot be run for cross-compile builds. Set FC flags to ensure 8-byte reals are being used.])])
AC_MSG_RESULT($test_result)

# If it failed, error out of configure with a helpful message.
if test $test_result = no; then
   AC_MSG_ERROR([Fortran compiler must be set to 8-byte reals. \
   Set FCFLAGS to allow this (-fdefault-real-8 for gfortran)])
fi

# Done testing Fortran compiler.
AC_LANG_POP(Fortran)
])

dnl Check that Fortran compiler is using 8-byte doubles.
dnl Ed Hartnett, 9/11/19
AC_DEFUN([GFDL_CHECK_FORTRAN_DOUBLE_8],
[
# Tell autoconf that language tests are in fortran.
AC_LANG_PUSH(Fortran)

# Check a fortran compile with a test program with this file
# extension.
AC_FC_SRCEXT(F90)

# Check the size of double.
AC_MSG_CHECKING([whether Fortran compiler has 8-byte double])
AC_RUN_IFELSE([AC_LANG_SOURCE([
        program foo
        double precision var
        if (sizeof(var) .ne. 8) then
        stop 2
        end if
        end program foo
        ])],
        [test_result=yes],
        [test_result=no],
        [AC_MSG_WARN([Test for 8-byte double cannot be run for cross-compile builds. Set FC flags to ensure 8-byte doubles are being used.])])
AC_MSG_RESULT($test_result)

# If it failed, error out of configure with a helpful message.
if test $test_result = no; then
   AC_MSG_ERROR([Fortran compiler must be set to 8-byte doubles. \
   Set FCFLAGS to allow this (-fdefault-double-8 for gfortran)])
fi

# Done testing Fortran compiler.
AC_LANG_POP(Fortran)
])


dnl Check that Fortran compiler can handle cray pointeres.
dnl Ed Hartnett, 9/11/19
AC_DEFUN([GFDL_CHECK_FORTRAN_CRAY_POINTERS],
[
# Tell autoconf that language tests are in fortran.
AC_LANG_PUSH(Fortran)

# Check a fortran compile with a test program with this file
# extension.
AC_FC_SRCEXT(F90)

# Check that code that uses a cray pointer can be compiled.
AC_MSG_CHECKING([whether Fortran compiler can handle cray pointers])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
        subroutine foo(bar)
        real :: dummy
        pointer( ptr, dummy )
        end subroutine foo
        ])],
        [test_result=yes],
        [test_result=no])
AC_MSG_RESULT($test_result)

# If it failed, error out of configure with a helpful message.
if test $test_result = no; then
   AC_MSG_ERROR([Fortran compiler must handle cray pointers. \
   Set FCFLAGS to allow this (-fcray-pointer for gfortran)])
fi

# Done testing Fortran compiler.
AC_LANG_POP(Fortran)
])
