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
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# The convenience library depends on its source.
noinst_LTLIBRARIES = libmemutils.la
libmemutils_la_SOURCES = memutils.F90

BUILT_SOURCES = memutils_mod.$(FC_MODEXT)
nodist_include_HEADERS = memutils_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
