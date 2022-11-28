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

# This is an automake file for the constants4 directory of the FMS
# package.

# Ryan Mulhall 7/22

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libconstants4.la

# The convenience library depends on its source.
libconstants4_la_SOURCES = \
  fmsconstantsr4.F90 \
  gfdl_constantsR4.h \
  gfs_constantsR4.h \
  geos_constantsR4.h \
  constantsr4.F90

fmsconstantsr4.$(FC_MODEXT): gfdl_constantsR4.h gfs_constantsR4.h geos_constantsR4.h
constantsr4_mod.$(FC_MODEXT): fmsconstantsr4.$(FC_MODEXT)

# Mod files are built and then installed as headers
MODFILES = \
  fmsconstantsr4.$(FC_MODEXT) \
  constantsr4_mod.$(FC_MODEXT)

BUILT_SOURCES = $(MODFILES)
nodist_include_HEADERS = $(MODFILES)

include $(top_srcdir)/mkmods.mk
