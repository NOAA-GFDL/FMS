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

# This is an automake file for the horiz_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/horiz_interp/include -I$(top_srcdir)/grid_utils
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libhoriz_interp.la

# Each convenience library depends on its source.
libhoriz_interp_la_SOURCES = \
  horiz_interp_bicubic.F90 \
  horiz_interp_bilinear.F90 \
  horiz_interp_conserve.F90 \
  horiz_interp.F90 \
  horiz_interp_spherical.F90 \
  horiz_interp_type.F90 \
  include/horiz_interp_bicubic.inc \
  include/horiz_interp_bilinear.inc \
  include/horiz_interp_conserve.inc \
  include/horiz_interp.inc \
  include/horiz_interp_spherical.inc \
  include/horiz_interp_type.inc \
  include/horiz_interp_bicubic_r4.fh \
  include/horiz_interp_conserve_xgrid.h \
  include/horiz_interp_conserve_xgrid.c \
  include/horiz_interp_bilinear_r4.fh \
  include/horiz_interp_conserve_r4.fh \
  include/horiz_interp_r4.fh \
  include/horiz_interp_spherical_r4.fh \
  include/horiz_interp_type_r4.fh \
  include/horiz_interp_bicubic_r8.fh \
  include/horiz_interp_bilinear_r8.fh \
  include/horiz_interp_conserve_r8.fh \
  include/horiz_interp_r8.fh \
  include/horiz_interp_spherical_r8.fh \
  include/horiz_interp_type_r8.fh

# Some mods are dependant on other mods in this dir.
horiz_interp_bicubic_mod.$(FC_MODEXT): horiz_interp_type_mod.$(FC_MODEXT)
horiz_interp_bilinear_mod.$(FC_MODEXT): horiz_interp_type_mod.$(FC_MODEXT)
horiz_interp_conserve_mod.$(FC_MODEXT): horiz_interp_type_mod.$(FC_MODEXT)
horiz_interp_spherical_mod.$(FC_MODEXT): horiz_interp_type_mod.$(FC_MODEXT)
horiz_interp_mod.$(FC_MODEXT): horiz_interp_bicubic_mod.$(FC_MODEXT) horiz_interp_type_mod.$(FC_MODEXT) \
                      horiz_interp_bilinear_mod.$(FC_MODEXT) horiz_interp_conserve_mod.$(FC_MODEXT) \
                      horiz_interp_spherical_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  horiz_interp_type_mod.$(FC_MODEXT) \
  horiz_interp_bicubic_mod.$(FC_MODEXT) \
  horiz_interp_bilinear_mod.$(FC_MODEXT) \
  horiz_interp_conserve_mod.$(FC_MODEXT) \
  horiz_interp_spherical_mod.$(FC_MODEXT) \
  horiz_interp_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
