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

# This is an automake file for the amip_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_interp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/fms2_io

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libamip_interp.la

# The convenience library depends on its source.
libamip_interp_la_SOURCES = amip_interp.F90

# Mod file depends on its o file, is built and then installed.
amip_interp_mod.mod: amip_interp.lo
BUILT_SOURCES = amip_interp_mod.mod
include_HEADERS = amip_interp_mod.mod

CLEANFILES = *.mod
