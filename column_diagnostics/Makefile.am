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
include/column_diagnostics.fh \
include/column_diagnostics_r4.inc \
include/column_diagnostics_r8.inc

column_diagnostics.$(FC_MOD_EXT):\
include/column_diagnostics.fh \
include/column_diagnostics_r4.inc \
include/column_diagnostics_r8.inc

BUILT_SOURCES = column_diagnostics_mod.$(FC_MODEXT)
nodist_include_HEADERS = column_diagnostics_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
