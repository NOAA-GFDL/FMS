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
#   Copyright (c) 2019 Seth Underwood <underwoo@underwoo.io>
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
#    -real_size 64: Intel compiler
#        -s real64: Cray
#              -r8: Portland Group compiler
#     -qrealsize=8: IBM compiler
AC_DEFUN([GX_FC_DEFAULT_REAL_KIND8_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran default REAL KIND 8 flag], [gx_cv_fc_default_real_kind8_flag],[
gx_cv_fc_default_real_kind8_flag=unknown
gx_fc_default_real_kind8_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fdefault-real-8' \
               '-real_size 64' \
               '-s real64' \
               '-r8' \
               '-qrealsize=8'; do
  test "x$ac_flag" != xnone && FCFLAGS="$gx_fc_default_real_kind8_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[
     program test
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
#    -real_size 32: Intel compiler
#        -s real32: Cray
#              -r4: Portland Group compiler
#     -qrealsize=4: IBM compiler
AC_DEFUN([GX_FC_DEFAULT_REAL_KIND4_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran default REAL KIND 4 flag], [gx_cv_fc_default_real_kind4_flag],[
gx_cv_fc_default_real_kind4_flag=unknown
gx_fc_default_real_kind4_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fdefault-real-4' \
               '-real_size 32' \
               '-s real32' \
               '-r4' \
               '-qrealsize=4'; do
  test "x$ac_flag" != xnone && FCFLAGS="$gx_fc_default_real_kind4_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[
     program test
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
AC_COMPILE_IFELSE([[
   program test
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
AC_CACHE_CHECK([for Fortran Cray Pointer flag], [gx_cv_fc_cray_ptr_flag],[
gx_cv_fc_cray_ptr_flag=unknown
gx_cray_ptr_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fcray-pointer' \
               '-Mcray=pointer'; do
  test "x$ac_flag" != xnone && FCFLAGS="$gx_cray_ptr_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[
     program test
     integer(kind=8) :: ipt
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
AC_CACHE_CHECK([if Fortran supports reading namelist from internal files], [gx_cv_fc_internal_file_nml],[
gx_cv_fc_internal_file_nml=unknown
AC_COMPILE_IFELSE([[
   program test
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
AC_LANG_PUSH([Fortran])
m4_ifval([$2],[gx_fc_check_mod_only=",only:$2"],[gx_fc_check_mod_only=""])
AS_LITERAL_WORD_IF([$1],
  [AS_VAR_PUSHDEF([gx_mod], [gx_cv_fc_check_mod_$1])],
  [AS_VAR_PUSHDEF([gx_mod], AS_TR_SH([gx_cv_check_mod_$1]))])
AC_CACHE_CHECK([for Fortran module $1], [gx_mod],
  [AC_COMPILE_IFELSE([[      program test
      use $1$gx_fc_check_mod_only
      end program test]],
    [AS_VAR_SET([gx_mod], [yes])],
    [AS_VAR_SET([gx_mod], [no])])])
AS_VAR_IF([gx_mod], [yes],
  [m4_default([$3],
    [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_MOD_$1), 1, [Define to 1 if the Fortran module $1 is found])])],
  [$4])
AS_VAR_POPDEF([gx_mod])
AC_LANG_POP([Fortran])
])
