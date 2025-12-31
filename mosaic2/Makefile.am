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

# This is an automake file for the mosaic directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/mosaic2/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.

noinst_LTLIBRARIES = libmosaic2.la

libmosaic2_la_SOURCES = \
mosaic2.F90 \
grid2.F90 \
include/mosaic2_r4.fh  include/mosaic2_r8.fh  include/mosaic2.inc \
include/grid2_r4.fh  include/grid2_r8.fh include/grid2.inc

# Some mods are dependant on other mods in this dir.
grid2_mod.$(FC_MODEXT): mosaic2_mod.$(FC_MODEXT) include/grid2_r4.fh include/grid2_r8.fh include/grid2.inc
mosaic2_mod.$(FC_MODEXT): include/mosaic2_r4.fh include/mosaic2_r8.fh include/mosaic2.inc

MODFILES = \
	mosaic2_mod.$(FC_MODEXT) \
	grid2_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
