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

# This is an automake file for the constants directory of the FMS
# package.

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/string_utils/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libstring_utils.la

# The convenience library depends on its source.
libstring_utils_la_SOURCES = \
  fms_string_utils.F90 \
  include/fms_string_utils.inc \
  include/fms_string_utils_r4.fh \
  include/fms_string_utils_r8.fh \
  fms_string_utils_binding.c

MODFILES = \
  fms_string_utils_mod.$(FC_MODEXT)
BUILT_SOURCES = $(MODFILES)
nodist_include_HEADERS = $(MODFILES)

include $(top_srcdir)/mkmods.mk
