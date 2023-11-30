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
#   GX_FC_ALLOW_ARG_MISMATCH([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
#
# DESCRIPTION
#
#   Set of functions that check for compiler flags to support legacy Fortran
#   language syntax
#
# LICENSE
#
#   Copyright (c) 2022 Seth Underwood <underwoo@underwoo.io>
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

# GX_FC_ALLOW_ARG_MISMATCH([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------
# Check if a compiler flag is required when calls to external procedures
# have argument mismatches between the different calls.
#
# Sets the variable FC_ALLOW_ARG_MISMATCH_FLAG to hold the flag.
#
# The known flags are:
# -fallow-argument-mismatch: gfortran (version 10 and later.  Not required
#                                      for versions less than 10.)
AC_DEFUN([GX_FC_ALLOW_ARG_MISMATCH],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran flag to allow procedure arg mismatch], [gx_cv_fc_allow_arg_mismatch_flag],[
gx_cv_fc_allow_arg_mismatch_flag=unknown
gx_fc_allow_arg_mismatch_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fallow-argument-mismatch'; do
test "x$ac_flag" != xnone && FCFLAGS="$gx_fc_allow_arg_mismatch_flag_FCFLAGS_save ${ac_flag}"
AC_COMPILE_IFELSE([[      program test
      logical(kind=8) :: arg8
      logical(kind=4) :: arg4
      call something(arg8)
      call something(arg4)
      end program test]],
      [gx_cv_fc_allow_arg_mismatch_flag=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$gx_fc_allow_arg_mismatch_flag_FCFLAGS_save
])
if test "x$gx_cv_fc_allow_arg_mismatch_flag" = xunknown; then
  m4_default([$2],
    [AC_MSG_ERROR([Unable to determine flag to allow argument mismatch])])
else
  FC_ALLOW_ARG_MISMATCH_FLAG=$gx_cv_fc_allow_arg_mismatch_flag
  if test "x$FC_ALLOW_ARG_MISMATCH_FLAG" = "xnone"; then
    FC_ALLOW_ARG_MISMATCH_FLAG=
  fi
  $1
fi
AC_LANG_POP([Fortran])
AC_SUBST([FC_ALLOW_ARG_MISMATCH_FLAG])
])
