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

# This is an automake file for the topography directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/topography/include

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
