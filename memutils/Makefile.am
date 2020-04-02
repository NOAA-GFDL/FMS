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

# This is an automake file for the memutils directory of the FMS package.

# Copyright 2019-2020, GFDL, Ed Hartnett, Seth Underwood

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

# The convenience library depends on its source.
libmemutils_la_SOURCES = memutils.F90
noinst_LTLIBRARIES = libmemutils.la

# Mod file depends on its o file, is built and then installed.
memutils_mod.mod: memutils.lo
BUILT_SOURCES = memutils_mod.mod
include_HEADERS = memutils_mod.mod

CLEANFILES = *.mod
