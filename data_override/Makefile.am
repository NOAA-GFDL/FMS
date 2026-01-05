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

# This is an automake file for the data_override directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/data_override/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libdata_override.la

# The convenience library depends on its source.
libdata_override_la_SOURCES = \
  get_grid_version.F90 \
  include/get_grid_version_r4.fh \
  include/get_grid_version_r8.fh \
  include/get_grid_version.inc \
  data_override.F90 \
  include/data_override_r4.fh \
  include/data_override_r8.fh \
  include/data_override.inc

# Some mods are dependent on other mods in this dir.
data_override_mod.$(FC_MODEXT): get_grid_version_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  data_override_mod.$(FC_MODEXT) \
  get_grid_version_mod.$(FC_MODEXT)

nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
