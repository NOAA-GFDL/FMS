#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
#***********************************************************************

# This is an automake file for the coupler directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/coupler/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcoupler.la

# Each convenience library depends on its source.
libcoupler_la_SOURCES = \
  coupler_types.F90 \
  ensemble_manager.F90 \
  atmos_ocean_fluxes.F90 \
  gex.F90 \
  include/coupler_types.inc \
  include/coupler_types_r4.fh \
  include/coupler_types_r8.fh

# Some mods are dependant on other mods in this dir.
atmos_ocean_fluxes_mod.$(FC_MODEXT): coupler_types_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  coupler_types_mod.$(FC_MODEXT) \
  ensemble_manager_mod.$(FC_MODEXT) \
  atmos_ocean_fluxes_mod.$(FC_MODEXT) \
  gex_mod.$(FC_MODEXT)

nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
