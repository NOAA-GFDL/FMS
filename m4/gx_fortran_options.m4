#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# ===========================================================================
#
# SYNOPSIS
#
#   GX_FC_DEFAULT_REAL_KIND8_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
#   GX_FC_DEFAULT_REAL_KIND4_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
#   GX_FC_QUAD_PRECISION()
#   GX_FC_CRAY_POINTER_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
#   GX_FC_INTERNAL_FILE_NML()
#   GX_FC_CHECK_MOD(module-name, [only], [action-if-found], [action-if-not-found])
#
# DESCRIPTION
#
#   Set of functions that check if the Fortran compiler supports certain
#   Fortran feature, or if a specific compiler flag is needed to support
#   the feature.  Full descriptions are avalable below.
#
# LICENSE
#
#   Copyright (c) 2019,2020 Seth Underwood <underwoo@underwoo.io>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.

# GX_FC_DEFAULT_REAL_KIND8_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------
# Look for the compiler flag that sets the default REAL kind to KIND=8.
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile with default REAL(KIND=8)) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.
#
# Sets the variable FC_DEFAULT_REAL_KIND8_FLAG to hold the flag.
#
# The known flags are:
# -fdefault-real-8: gfortran
#    -real-size 64: Intel compiler
#    -real_size 64: Intel compiler
#        -s real64: Cray
#              -r8: Portland Group compiler
#     -qrealsize=8: IBM compiler
AC_DEFUN([GX_FC_DEFAULT_REAL_KIND8_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran flag needed to accept default REAL(KIND=8)], [gx_cv_fc_default_real_kind8_flag],[
gx_cv_fc_default_real_kind8_flag=unknown
gx_fc_default_real_kind8_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fdefault-real-8' \
               '-real-size 64' \
               '-real_size 64' \
               '-s real64' \
               '-r8' \
               '-qrealsize=8'; do
  test "x$ac_flag" != xnone && FCFLAGS="$gx_fc_default_real_kind8_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[      program test
      interface
      subroutine test_sub(a)
      real(kind=selected_real_kind(15,307)) :: a
      end subroutine test_sub
      end interface
      real :: b=1.0
      call test_sub(b)
      end program test]],
     [gx_cv_fc_default_real_kind8_flag=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$gx_fc_default_real_kind8_flag_FCFLAGS_save
])
if test "x$gx_cv_fc_default_real_kind8_flag" = xunknown; then
  m4_default([$2],
              [AC_MSG_ERROR([Fortran cannot set default real kind to 8])])
else
  FC_DEFAULT_REAL_KIND8_FLAG=$gx_cv_fc_default_real_kind8_flag
  if test "x$FC_DEFAULT_REAL_KIND8_FLAG" = xnone; then
    FC_DEFAULT_REAL_KIND8_FLAG=
  fi
  $1
fi
AC_LANG_POP([Fortran])
AC_SUBST([FC_DEFAULT_REAL_KIND8_FLAG])
])

# GX_FC_DEFAULT_REAL_KIND4_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------
# Look for the compiler flag that sets the default REAL kind to KIND=4.
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile with default REAL(KIND=4)) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.
#
# Sets the variable FC_DEFAULT_REAL_KIND4_FLAG to hold the flag.
#
# The known flags are:
#             none: gfortran (gfortran does not have an option to set the
#                   default REAL kind to KIND=4)
#    -real-size 32: Intel compiler
#    -real_size 32: Intel compiler
#        -s real32: Cray
#              -r4: Portland Group compiler
#     -qrealsize=4: IBM compiler
AC_DEFUN([GX_FC_DEFAULT_REAL_KIND4_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran flag needed to accept default REAL(KIND=4)], [gx_cv_fc_default_real_kind4_flag],[
gx_cv_fc_default_real_kind4_flag=unknown
gx_fc_default_real_kind4_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fdefault-real-4' \
               '-real-size 32' \
               '-real_size 32' \
               '-s real32' \
               '-r4' \
               '-qrealsize=4'; do
  test "x$ac_flag" != xnone && FCFLAGS="$gx_fc_default_real_kind4_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[      program test
      interface
      subroutine test_sub(a)
      real(kind=selected_real_kind(6, 37)) :: a
      end subroutine test_sub
      end interface
      real :: b=1.0
      call test_sub(b)
      end program test]],
     [gx_cv_fc_default_real_kind4_flag=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$gx_fc_default_real_kind4_flag_FCFLAGS_save
])
if test "x$gx_cv_fc_default_real_kind4_flag" = xunknown; then
  m4_default([$2],
              [AC_MSG_ERROR([Fortran cannot set default real kind to 4])])
else
  FC_DEFAULT_REAL_KIND4_FLAG=$gx_cv_fc_real_kind4_flag
  if test "x$FC_DEFAULT_REAL_KIND4_FLAG" = xnone; then
    FC_DEFAULT_REAL_KIND4_FLAG=
  fi
  $1
fi
AC_LANG_POP([Fortran])
AC_SUBST([FC_DEFAULT_REAL_KIND4_FLAG])
])

# GX_FC_QUAD_PRECISION
# -----------------------------------------------------------------------------
# Determine if the Fortran compiler and target system have support for IEEE 754,
# quadruple precision.  If supported, sets the define HAVE_QUAD_PRECISION.
AC_DEFUN([GX_FC_QUAD_PRECISION],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([if Fortran and target have IEEE 754 support], [gx_cv_fc_quad_precision],[
gx_cv_fc_quad_precision=unknown
AC_COMPILE_IFELSE([[      program test
      real(KIND=selected_real_kind(33, 4931)) :: quad
      end program test]],
   [gx_cv_fc_quad_precision=yes],
   [gx_cv_fc_quad_precision=no])])
if test "x$gx_cv_fc_quad_precision" = "xyes"; then
   AC_DEFINE([HAVE_QUAD_PRECISION], 1,
             [Define to 1 if your Fortran and system have IEEE 754 support])
fi
AC_LANG_POP([Fortran])
])

# GX_FC_CRAY_POINTER_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Look for the compiler flag that allows Fortran Cray Pointers.  Cray
# pointers are an are part of a non-standard extension that provides a
# C-like pointer in Fortran.  Call ACTION-IF-SUCCESS (defaults to
# nothing) if successful (i.e. can use Cray pointers) and
# ACTION-IF-FAILURE (defaults to failing with an error message) if not.
#
# Sets the variable FC_CRAY_POINTER_FLAG to hold the flag, and defines
# HAVE_CRAY_POINTER.
#
# The known flags are:
# -fcray-pointer: gfortran
#           none: Intel compiler (No option required for Cray Pointers)
#        unknown: Cray
# -Mcray=pointer: Portland Group compiler
#           none: IBM compiler (No option required for Cray Pointers)
AC_DEFUN([GX_FC_CRAY_POINTER_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran flag needed to accept Cray pointers], [gx_cv_fc_cray_ptr_flag],[
gx_cv_fc_cray_ptr_flag=unknown
gx_cray_ptr_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fcray-pointer' \
               '-Mcray=pointer'; do
  test "x$ac_flag" != xnone && FCFLAGS="$gx_cray_ptr_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[      program test
      integer iarri(10)
      pointer (ipt, iarr)
      end program test]],
     [gx_cv_fc_cray_ptr_flag=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$gx_cray_ptr_flag_FCFLAGS_save
])
if test "x$gx_cv_fc_cray_ptr_flag" = "xunknown"; then
  m4_default([$2],
              [AC_MSG_ERROR([Fortran cannot use Cray pointers])])
else
  AC_DEFINE([HAVE_CRAY_POINTER], 1,
            [Define to 1 if your Fortran compiler supports cray pointers])
  FC_CRAY_POINTER_FLAG=$gx_cv_fc_cray_ptr_flag
  if test "x$FC_CRAY_POINTER_FLAG" = xnone; then
    FC_CRAY_POINTER_FLAG=
  fi
  $1
fi
AC_LANG_POP([Fortran])
AC_SUBST([FC_CRAY_POINTER_FLAG])
])

# GX_FC_INTERNAL_FILE_NML
# -----------------------------------------------------------------------------
# Determine if the Fortran compiler supports reading Fortran namelists from
# an internal file.  If supported, sets the define HAVE_INTERNAL_NML.
AC_DEFUN([GX_FC_INTERNAL_FILE_NML],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([if $[]_AC_FC[] supports reading namelists from internal files], [gx_cv_fc_internal_file_nml],[
gx_cv_fc_internal_file_nml=unknown
AC_COMPILE_IFELSE([[      program test
      implicit none
      integer :: a = 1
      real :: b = 0.1
      character(LEN=20) :: internal_nml ="&test_nml a=2 b=1.0/"
      namelist /test_nml/ a, b
      read(internal_nml,test_nml)
      end program test]],
   [gx_cv_fc_internal_file_nml=yes],
   [gx_cv_fc_internal_file_nml=no])])
if test "x$gx_cv_fc_internal_file_nml" = "xyes"; then
   AC_DEFINE([HAVE_INTERNAL_NML], 1,
             [Define to 1 if your Fortran compiler supports reading namelists
              from internal files])
fi
AC_LANG_POP([Fortran])
])

# GX_FC_CHECK_MOD(module-name, [only], [action-if-found], [action-if-not-found])
# -----------------------------------------------------------------------------
# Check if a Fortran module module-name is available.  Execute shell commands
# action-if-found, otherwise execute action-if-not-found.  If only is specified
# then check if the Fortran module has the given symbol.
AC_DEFUN([GX_FC_CHECK_MOD],[
_AC_FORTRAN_ASSERT()dnl
m4_ifval([$2],[gx_fc_check_mod_only=",only:$2"],[gx_fc_check_mod_only=""])
AS_LITERAL_WORD_IF([$1],
  [AS_VAR_PUSHDEF([gx_mod], [gx_cv_fc_check_mod_$1])],
  [AS_VAR_PUSHDEF([gx_mod], AS_TR_SH([gx_cv_check_mod_$1]))])
dnl Autoconf does not pass CPPFLAGS to the Fortran tests.  As the user may
dnl define the include options in CPPFLAGS instead of FFLAGS or FCFLAGS, we
dnl force CPPFLAGS te be part of FFLAGS/FCFLAGS
gx_save_[]_AC_LANG_PREFIX[]FLAGS=${[]_AC_LANG_PREFIX[]FLAGS}
_AC_LANG_PREFIX[]FLAGS="${[]_AC_LANG_PREFIX[]FLAGS} $CPPFLAGS"
AC_CACHE_CHECK([for Fortran module $1], [gx_mod],
  [AC_COMPILE_IFELSE([[      program test
      use $1$gx_fc_check_mod_only
      end program test]],
    [AS_VAR_SET([gx_mod], [yes])],
    [AS_VAR_SET([gx_mod], [no])])])
_AC_LANG_PREFIX[]FLAGS="${gx_save_[]_AC_LANG_PREFIX[]FLAGS}"
AS_VAR_IF([gx_mod], [yes],
  [m4_default([$3],
    [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_MOD_$1), 1, [Define to 1 if the Fortran module $1 is found])])],
  [$4])
AS_VAR_POPDEF([gx_mod])
])#GX_FC_CHECK_MOD

# GX_FORTRAN_CHECK_HEADERS(header, [action-if-found], [action-if-not-found])
# -----------------------------------------------------------------------
# Check if a Fortran include file is available.
AC_DEFUN([GX_FORTRAN_CHECK_HEADERS], [
_AC_FORTRAN_ASSERT()dnl
AS_LITERAL_WORD_IF([$1],
  [AS_VAR_PUSHDEF([gx_header], [gx_cv_fortran_check_headers_$1])],
  [AS_VAR_PUSHDEF([gx_header], [AS_TR_SH([gx_cv_fortran_check_headers_$1])])])
dnl Ensure the Fortran compiler will run the preprocessor
_GX_FORTRAN_PP_SRCEXT_PUSH([F])
dnl Autoconf does not pass CPPFLAGS to the Fortran tests.  As the user may
dnl define the include options in CPPFLAGS instead of FFLAGS or FCFLAGS, we
dnl force CPPFLAGS te be part of FFLAGS/FCFLAGS
gx_save_[]_AC_LANG_PREFIX[]FLAGS=${[]_AC_LANG_PREFIX[]FLAGS}
_AC_LANG_PREFIX[]FLAGS="${[]_AC_LANG_PREFIX[]FLAGS} $CPPFLAGS"
AC_CACHE_CHECK([for $1 usability], [gx_header],
  [AC_COMPILE_IFELSE(AC_LANG_PROGRAM([],
    [@%:@include <$1>]),
    [AS_VAR_SET([gx_header], [yes])],
    [AS_VAR_SET([gx_header], [no])])])
_AC_LANG_PREFIX[]FLAGS="${gx_save_[]_AC_LANG_PREFIX[]FLAGS}"
_GX_FORTRAN_PP_SRCEXT_POP([F])
AS_VAR_IF([gx_header], [yes],
  [m4_default([$2],
    [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_$1), 1, [Define to 1 if the Fortran include file $1 is found])])],
  [$3])
AS_VAR_POPDEF([gx_header])
])# GX_FC_CHECK_HEADERS

# GX_FORTRAN_SEARCH_LIBS(FUNCTION, SEARCH-LIBS, [CALL_PREFIX], [CALL_SYNTAX],
#                        [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#                        [OTHER-LIBRARIES])
# ------------------------------------------------------------------------------
# Search for a Fortran library defining FUNCTION, if it's not already availabe.
#
# This expands AC_SEARCH_LIBS for Fortran as the AC function does not work in
# all cases for Fortran.
AC_DEFUN([GX_FORTRAN_SEARCH_LIBS], [
_AC_FORTRAN_ASSERT()dnl
AS_VAR_PUSHDEF([gx_search], [gx_cv_fortran_search_libs_$1])
dnl Ensure the Fortran compiler will run the preprocessor
_GX_FORTRAN_PP_SRCEXT_PUSH([F])
dnl Autoconf does not pass CPPFLAGS to the Fortran tests.  As the user may
dnl define the include options in CPPFLAGS instead of FFLAGS or FCFLAGS, we
dnl force CPPFLAGS te be part of FFLAGS/FCFLAGS
gx_save_FLAGS="${[]_AC_LANG_PREFIX[]FLAGS}"
_AC_LANG_PREFIX[]FLAGS="${[]_AC_LANG_PREFIX[]FLAGS} $CPPFLAGS"
dnl Prepare the prefix and syntax sections so they will compile correctly
dnl with the Fortran compiler
m4_foreach(gx_line, m4_split($3, m4_newline()),
  [m4_append([_gx_fortran_sanitized_call_prefix], m4_bmatch(gx_line, [^\w.*], [      ]gx_line, [gx_line]), m4_newline())dnl
])
gx_fortran_sanitized_call_prefix="_gx_fortran_sanitized_call_prefix"
gx_fortran_func_line="m4_bregexp([$4], [^.*$], [      \&])"
AC_CACHE_CHECK([for library containing $1], [gx_search],
[gx_fortran_search_libs_save_LIBS=$LIBS
AC_LANG_CONFTEST([AC_LANG_PROGRAM([],[[$gx_fortran_sanitized_call_prefix
$gx_fortran_func_line]])])
dnl Search for the library
for gx_lib in '' $2; do
  if test -z "$gx_lib"; then
    gx_res="none required"
  else
    gx_res=-l$gx_lib
    LIBS="$gx_res $7 $gx_fortran_search_libs_save_LIBS"
  fi
  AC_LINK_IFELSE([], [AS_VAR_SET([gx_search], [$gx_res])])
  AS_VAR_SET_IF([gx_search], [break])
done
AS_VAR_SET_IF([gx_search], , [AS_VAR_SET([gx_search], [no])])
rm conftest.$ac_ext
LIBS=$gx_fortran_search_libs_save_LIBS])
AS_VAR_COPY([gx_res], [gx_search])
AS_IF([test "$gx_res" != no],
  [test "$gx_res" = "none required" || LIBS="$gx_res $LIBS"
  $5],
  [$6])
_AC_LANG_PREFIX[]FLAGS="$gx_save_FLAGS"
_GX_FORTRAN_PP_SRCEXT_POP([F])
AS_VAR_POPDEF([gx_search])
m4_ifdef([gx_line], [m4_undefine([gx_line])])
m4_ifdef([_gx_fortran_sanitized_call_prefix], [m4_undefine([_gx_fortran_sanitized_call_prefix])])
])# GX_FORTRAN_SEARCH_LIBS

# _GX_FORTRAN_PP_SRCEXT_PUSH(EXT) will ensure the Fortran test extension (stored
# in ac_ext) will cause the Fortran compiler to preprocess the test source file.
# Most Fortran compilers will preprocess the file based on the file extension,
# and of the known extension, F appears to work for all compilers.  This
# function accepts the file extension (without the preceeding .), similar
# to the AC_FC_SRCEXT and AC_FC_PP_SRCEXT macros.  If the extension is not
# provided, the default is 'F'.
# Unfortunately, this will not work for compilers that require a specific flag.
AC_DEFUN([_GX_FORTRAN_PP_SRCEXT_PUSH],[
_gx_fortran_pp_srcext_save=$ac_ext
AS_VAR_SET_IF([1], [ac_ext=$1], [ac_ext="F"])
])# _GX_FORTRAN_PP_SRCEXT_PUSH

# _GX_FORTRAN_PP_SRCEXT_POP() reverses the extension change done in
# _GX_FORTRAN_PP_SRCEXT_PUSH.  If this pop is called without a preceeding
# push, the default Fortran file extension (f) will be used.
AC_DEFUN([_GX_FORTRAN_PP_SRCEXT_POP], [
AS_VAR_SET_IF([_gx_fortran_pp_srcext_save], [ac_ext=${_gx_fortran_pp_srcext_save}], [ac_ext="f"])
unset _gx_fortran_pp_srcext_save
])# _GX_FORTRAN_PP_SRCEXT_POP
