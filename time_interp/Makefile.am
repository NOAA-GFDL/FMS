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

# This is an automake file for the time_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/time_interp/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtime_interp.la

# Each convenience library depends on its source.
libtime_interp_la_SOURCES = \
  time_interp.F90 \
  time_interp_external2.F90 \
  include/time_interp_r4.fh \
  include/time_interp_r8.fh \
  include/time_interp.inc \
  include/time_interp_external2_bridge_r4.fh \
  include/time_interp_external2_bridge_r8.fh \
  include/time_interp_external2_bridge.inc \
  include/time_interp_external2_r4.fh \
  include/time_interp_external2_r8.fh \
  include/time_interp_external2.inc

# Some mods are dependant on other mods in this dir.
time_interp_external2_mod.$(FC_MODEXT): time_interp_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  time_interp_mod.$(FC_MODEXT) \
  time_interp_external2_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
