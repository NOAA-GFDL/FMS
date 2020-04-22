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
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/data_override
AM_CPPFLAGS += -I${top_builddir}/field_manager
AM_CPPFLAGS += -I${top_builddir}/platform

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcoupler_types.la libensemble_manager.la \
libatmos_ocean_fluxes.la

# Each convenience library depends on its source.
libcoupler_types_la_SOURCES = coupler_types.F90
libensemble_manager_la_SOURCES = ensemble_manager.F90
libatmos_ocean_fluxes_la_SOURCES = atmos_ocean_fluxes.F90

# Each mod file depends on the .o file.
coupler_types_mod.mod: coupler_types.lo
ensemble_manager_mod.mod: ensemble_manager.lo
atmos_ocean_fluxes_mod.mod: atmos_ocean_fluxes.lo

# Some mods are dependant on other mods in this dir.
atmos_ocean_fluxes.lo: coupler_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = coupler_types_mod.mod ensemble_manager_mod.mod \
atmos_ocean_fluxes_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
