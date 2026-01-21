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
#   GX_FC_CLASS_CHAR_ARRAY_BUG_CHECK([ACTION-IF-BUG-PRESENT = FAILURE])
#
# DESCRIPTION
#
#   Set of functions to check if the compilers have any known bugs.
#   Full descriptions are available below.
#
# LICENSE
#
#   Copyright (c) 2020 Seth Underwood <underwoo@underwoo.io>, @uramirez8707
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

# GX_FC_CLASS_CHAR_ARRAY_BUG_CHECK([ACTION-IF-BUG-PRESENT = FAILURE])
# ----------------------------------------------------------------------
# Check if the Fortran compiler supports has the class character array
# assign bug.  If the Fortran compiler does have the bug, call
# ACTION-IF-BUG-PRESENT (defaults to failing).
AC_DEFUN([GX_FC_CLASS_CHAR_ARRAY_BUG_CHECK],[
_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([if $[]_AC_FC[] has the class character array assign bug], [gx_cv_class_char_array_bug_check],[dnl
gx_cv_class_char_array_bug_check=yes
AC_COMPILE_IFELSE([[      subroutine test_sub(ctype)
        class(*), intent(out) :: ctype

        select type(ctype)
          type is (character(len=*))
            ctype(:) = ""
        end select
      end subroutine test_sub]],
     [gx_cv_class_char_array_bug_check=no])
])
AS_IF([test "x$gx_cv_class_char_array_bug_check" = xyes],[dnl
  AC_DEFINE([HAVE_CLASS_CHAR_ARRAY_BUG], 1,
    [Define to 1 if the Fortran compiler has the class character array bug])
  m4_default([$1],
    [AC_MSG_ERROR([The Fortran compiler has the class, character array assing bug.  libFMS cannot be built with this compiler.])])
])dnl
])
