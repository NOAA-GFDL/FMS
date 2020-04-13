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

# This is an automake file for the time_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms2_io
AM_CPPFLAGS += -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/horiz_interp

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtime_interp.la libtime_interp_external.la	\
libtime_interp_external2.la

# Each convenience library depends on its source.
libtime_interp_la_SOURCES = time_interp.F90
libtime_interp_external_la_SOURCES = time_interp_external.F90
libtime_interp_external2_la_SOURCES = time_interp_external2.F90

# Each mod file depends on the .o file.
time_interp_mod.mod: time_interp.lo
time_interp_external_mod.mod: time_interp_external.lo
time_interp_external2_mod.mod: time_interp_external2.lo

# Some mods are dependant on other mods in this dir.
time_interp_external.lo: time_interp_mod.mod
time_interp_external2.lo: time_interp_mod.mod

# Mod files are built and then installed as headers.
MODFILES = time_interp_mod.mod time_interp_external_mod.mod	\
time_interp_external2_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
