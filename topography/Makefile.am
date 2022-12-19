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
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/topography -I$(top_srcdir)/topography/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtopography.la

# Each convenience library depends on its source.
libtopography_la_SOURCES = \
    topography.F90 \
    include/topography_r4.fh \
    include/topography_r8.fh \
    include/topography.inc \
    gaussian_topog.F90 \
    include/gaussian_topog_r4.fh \
    include/gaussian_topog_r8.fh \
    include/gaussian_topog.inc

gaussian_topog_mod.$(FC_MODEXT): \
    include/gaussian_topog_r4.fh \
    include/gaussian_topog_r8.fh \
    include/gaussian_topog.inc

topography_mod.$(FC_MODEXT): \
    gaussian_topog_mod.$(FC_MODEXT) \
    include/topography_r4.fh \
    include/topography_r8.fh \
    include/topography.inc

# Mod files are built and then installed as headers.
MODFILES = \
  gaussian_topog_mod.$(FC_MODEXT) \
  topography_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
