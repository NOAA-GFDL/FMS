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

# This is an automake file for the field_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/fms2_io

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfield_manager.la libfm_util.la

# Each convenience library depends on its source.
libfield_manager_la_SOURCES = field_manager.F90
libfm_util_la_SOURCES = fm_util.F90

# Each mod file depends on the .o file.
field_manager_mod.mod: field_manager.lo
fm_util_mod.mod: fm_util.lo

# Some mods are dependant on other mods in this dir.
fm_util.lo: field_manager_mod.mod

# Mod files are built and then installed as headers.
MODFILES = field_manager_mod.mod fm_util_mod.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

EXTRA_DIST = parse.inc

CLEANFILES = *.mod
