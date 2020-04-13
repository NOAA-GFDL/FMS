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

# This is an automake file for the random_numbers directory of the FMS
# package.

# Ed Hartnett 2/28/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_manager

# Build these uninstalled convenience library.
noinst_LTLIBRARIES = libMersenneTwister.la librandom_numbers.la

# Each convenience library depends on its source.
librandom_numbers_la_SOURCES = random_numbers.F90
libMersenneTwister_la_SOURCES = MersenneTwister.F90

# Each mod file depends on the .o file.
random_numbers_mod.mod: random_numbers.lo
mersennetwister_mod.mod: MersenneTwister.lo

# Some mods are dependant on other mods in this dir.
random_numbers.lo: mersennetwister_mod.mod

# Mod files are built and then installed as headers.
MODFILES = mersennetwister_mod.mod random_numbers_mod.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

CLEANFILES = *.mod
