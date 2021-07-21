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
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.

noinst_LTLIBRARIES = libmosaic.la

libmosaic_la_SOURCES = \
mosaic.F90 \
grid.F90 \
gradient.F90 \
create_xgrid.c \
gradient_c2l.c \
interp.c \
mosaic_util.c \
read_mosaic.c \
constant.h \
create_xgrid.h \
gradient_c2l.h \
interp.h \
mosaic_util.h \
read_mosaic.h

# Some mods are dependant on other mods in this dir.
grid_mod.$(FC_MODEXT): mosaic_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  mosaic_mod.$(FC_MODEXT) \
  grid_mod.$(FC_MODEXT) \
  gradient_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
