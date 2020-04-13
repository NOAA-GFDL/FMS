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

# This is an automake file for the mosaic directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/fms2_io

# Build these uninstalled convenience libraries.

noinst_LTLIBRARIES = libmosaic.la libmosaic2.la libgrid.la libgradient.la \
libmosaic_c.la

# Each convenience library depends on its source. Lump all the C code
# into libmosaic_c.la.
libmosaic_la_SOURCES = mosaic.F90
libmosaic2_la_SOURCES = mosaic2.F90
libgrid_la_SOURCES = grid.F90
libgradient_la_SOURCES = gradient.F90
libmosaic_c_la_SOURCES = create_xgrid.c gradient_c2l.c interp.c		\
mosaic_util.c read_mosaic.c constant.h create_xgrid.h gradient_c2l.h	\
interp.h mosaic_util.h read_mosaic.h

# Each mod file depends on the .o file.
mosaic_mod.mod: mosaic.lo
mosaic2_mod.mod: mosaic2.lo
grid_mod.mod: grid.lo
gradient_mod.mod: gradient.lo

# Some mods are dependant on other mods in this dir.
grid.lo: mosaic_mod.mod

# Mod files are built and then installed as headers.
MODFILES = mosaic_mod.mod mosaic2_mod.mod grid_mod.mod gradient_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
