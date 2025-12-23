#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
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
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#   Copyright (c) 2022 Seth Underwood <underwoo@underwoo.io>

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
