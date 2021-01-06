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

# This is an automake file for the coupler directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcoupler.la

# Each convenience library depends on its source.
libcoupler_la_SOURCES = \
  coupler_types.F90 \
  ensemble_manager.F90 \
  atmos_ocean_fluxes.F90

# Some mods are dependant on other mods in this dir.
atmos_ocean_fluxes_mod.$(FC_MODEXT): coupler_types_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  coupler_types_mod.$(FC_MODEXT) \
  ensemble_manager_mod.$(FC_MODEXT) \
  atmos_ocean_fluxes_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
