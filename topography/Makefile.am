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

# This is an automake file for the topography directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/fms2_io

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libgaussian_topog.la libtopography.la

# Each convenience library depends on its source.
libtopography_la_SOURCES = topography.F90
libgaussian_topog_la_SOURCES = gaussian_topog.F90

# Each mod file depends on the .o file.
gaussian_topog_mod.mod: gaussian_topog.lo
topography_mod.mod: topography.lo

# Some mods are dependant on other mods in this dir.
topography.lo: gaussian_topog_mod.mod

# Mod files are built and then installed as headers.
MODFILES = gaussian_topog_mod.mod topography_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
