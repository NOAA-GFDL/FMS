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
#   GX_MPI([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
#   GX_MPI_FC_LEGACY_INTERFACE()
#
# DESCRIPTION
#
#   Determine the compiler can find the MPI headers and library.
#   This required gx_fortran_options.m4 for GX_FC_CHECK_MOD.
#
#   Also test if the MPI library uses legacy, Fortran interfaces.  In some compilers
#   these legacy interfaces can lead to an error with non-matching arguments.
#   In particular, GCC version >= 11.0.0
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

# GX_MPI([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------
# Check if the compiler can find the MPI headers and libraries
AC_DEFUN([GX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AS_VAR_PUSHDEF([gx_cv_mpi_h],[gx_cv_[]_AC_LANG_ABBREV[]_mpi_h])
dnl Check for the MPI header.  Here we use AC_TRY_COMPILE as AC_CHECK_HEADER
dnl Note, we use AC_TRY_COMPILE as AC_CHECK_HEADER will call $CPP. Since
dnl CC may be a mpi-specific compiler (e.g. mpicc), we don't want to use $CPP.
AC_LANG_CASE([C], [
    AC_REQUIRE([AC_PROG_CC])
    AC_CACHE_CHECK([for mpi.h], [gx_cv_mpi_h], [dnl
    gx_cv_mpi_h=no
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <mpi.h>])], [dnl
      gx_cv_mpi_h=yes])])
],
[C++], [
    AC_REQUIRE([AC_PROG_CXX])
    AC_CACHE_CHECK([for mpi.h], [gx_cv_mpi_h], [dnl
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <mpi.h>])],[dnl
      gx_cv_mpi_h=yes])])
],
[Fortran 77], [
    AC_REQUIRE([AC_PROG_F77])
    AC_CACHE_CHECK([for mpif.h], [gx_cv_mpi_h], [dnl
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [      include "mpif.h"])], [dnl
      gx_cv_mpi_h=yes])])
],
[Fortran], [
    AC_REQUIRE([AC_PROG_FC])
    AC_CACHE_CHECK([for mpif.h], [gx_cv_mpi_h], [dnl
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [      include "mpif.h"])], [dnl
      gx_cv_mpi_h=yes])])
    GX_FC_CHECK_MOD([mpi])
])
# Check for library
AC_SEARCH_LIBS([MPI_Init], [mpi])

AS_VAR_PUSHDEF([gx_mpi_lang_usable], [gx_[]_AC_LANG_ABBREV[]mpi])
AC_LANG_CASE([C], [test "$gx_cv_mpi_h" = yes -a "x$ac_cv_search_MPI_Init" != "xno" && gx_mpi_lang_usable=yes],
    [C++], [test "$gx_cv_mpi_h" = yes -a "x$ac_cv_search_MPI_Init" != "xno" && gx_mpi_lang_usable=yes],
    [Fortran 77], [test "$gx_cv_mpi_h" = yes -a "x$ac_cv_search_MPI_Init" != "xno" && gx_mpi_lang_usable=yes],
    [Fortran], [test \( "$gx_cv_mpi_h" = yes -o $gx_cv_fc_check_mod_mpi = yes \) -a "x$ac_cv_search_MPI_Init" != "xno" && gx_mpi_lang_usable=yes
])
AS_VAR_IF([gx_mpi_lang_usable], [yes], [dnl
    m4_default([$1], [])
    AC_LANG_CASE([C], [AC_DEFINE([HAVE_MPI_H], 1, [Define to 1 if the mpi.h header file $1 is found])],
        [C++], [AC_DEFINE([HAVE_MPI_H], 1, [Define to 1 if the mpi.h header file $1 is found])],
        [Fortran 77], [AC_DEFINE([HAVE_MPIF_H], 1, [Define to 1 if the mpif.h header file $1 is found])],
        [Fortran], [AC_DEFINE([HAVE_MPIF_H], 1, [Define to 1 if the mpif.h header file $1 is found])]
    )], [dnl
    m4_default([$2], [AC_MSG_ERROR([Unable to find the MPI headers or library for _AC_LANG])])
])
AS_VAR_POPDEF([gx_mpi_lang_usable])
AS_VAR_POPDEF([gx_cv_mpi_h])
])

# GX_MPI_FC_LEGACY_INTERFACE()
# ----------------------------------------------------------------------
# Check if the Fortran MPI library uses legacy interfaces.
#
# If the Fortran MPI library uses legacy interfaces, the variable
# HAVE_MPI_FC_LEGACY will be set.
AC_DEFUN([GX_MPI_FC_LEGACY_INTERFACE], [
    AC_LANG_ASSERT([Fortran])
    AC_MSG_CHECKING([if the MPI Fortran library uses legacy interfaces])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [dnl
      use mpi
      integer(kind=8) :: rarg8(10), sarg8(10)
      integer(kind=4) :: rarg4(10), sarg4(10)
      integer :: ssize(10), rsize(10), sdispl(10), rdispl(10), ierr
      call MPI_Alltoallv(sarg4, ssize, sdispl, MPI_INTEGER4, &
            rarg4, rsize, rdispl, MPI_INTEGER4, 1, error)
      call MPI_Alltoallv(sarg8, ssize, sdispl, MPI_INTEGER8, &
            rarg8, rsize, rdispl, MPI_INTEGER8, 1, error)])], [dnl
            AC_MSG_RESULT(no)], [dnl
            HAVE_MPI_FC_LEGACY=yes
            AC_SUBST([HAVE_MPI_FC_LEGACY])
            AC_MSG_RESULT(yes)]
    )
])
