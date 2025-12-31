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

# This is an automake file for the column_diagnostics directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/column_diagnostics/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libcolumn_diagnostics.la

# The convenience library depends on its source.
libcolumn_diagnostics_la_SOURCES = column_diagnostics.F90 \
include/column_diagnostics.inc \
include/column_diagnostics_r4.fh \
include/column_diagnostics_r8.fh

column_diagnostics.$(FC_MOD_EXT):\
include/column_diagnostics.inc \
include/column_diagnostics_r4.fh \
include/column_diagnostics_r8.fh

BUILT_SOURCES = column_diagnostics_mod.$(FC_MODEXT)
nodist_include_HEADERS = column_diagnostics_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
