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
#   GX_FC_08_OPEN_NEWUNIT([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
#
# DESCRIPTION
#
#   Set of functions that check if the Fortran compiler supports certain
#   Fortran 2008 features, or if a specific compiler flag is needed to support
#   the feature.  Full descriptions are avalable below.
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
#   Copyright (c) 2020 Seth Underwood <underwoo@underwoo.io>

# GX_FC_08_OPEN_NEWUNIT([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# ----------------------------------------------------------------------
# Check if the Fortran compiler supports OPEN(NEWUNIT=iunit,...).  The
# NEWUNIT option for OPEN will select, automatically, a new free unit
# number.  If the compiler supports this, the preprocessor macro
# HAVE_OPEN_NEWUNIT will be set.
AC_DEFUN([GX_FC_08_OPEN_NEWUNIT],[
_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([if $[]_AC_FC[] supports OPEN(NEWUNIT=iunit,...)], [gx_cv_fc_08_open_newfile],[dnl
gx_cv_fc_08_open_newfile=no
AC_COMPILE_IFELSE([[      program test
      integer :: funit
      open(NEWUNIT=funit, STATUS='SCRATCH')
      close(funit)
      end program test]],
     [gx_cv_fc_08_open_newfile=yes])
])
AS_IF([test "x$gx_cv_fc_08_open_newfile" = xno],[$2],[dnl
  AC_DEFINE([HAVE_OPEN_NEWUNIT], 1,
            [Define to 1 if the Fotran compiler supports OPEN(NEWUNIT=iunit,...)])
  $1
])dnl
])# GC_FC_08_OPEN_NEWUNIT
