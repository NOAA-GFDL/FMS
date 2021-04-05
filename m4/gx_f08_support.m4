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
#   Copyright (c) 2020 Seth Underwood <underwoo@underwoo.io>
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
