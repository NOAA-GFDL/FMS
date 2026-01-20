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

# This is an automake file for the tridiagonal directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/tridiagonal/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libtridiagonal.la

# The convenience library depends on its source.
libtridiagonal_la_SOURCES = tridiagonal.F90           \
                            include/tridiagonal.inc   \
                            include/tridiagonal_r4.fh \
                            include/tridiagonal_r8.fh

# Mod file depends on its o file, is built and then installed.
tridiagonal.lo: tridiagonal_mod.$(FC_MODEXT)

BUILT_SOURCES = tridiagonal_mod.$(FC_MODEXT)
nodist_include_HEADERS = tridiagonal_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
