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

# This is an automake file for the astronomy directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/astronomy/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libastronomy.la

# The convenience library depends on its source.
libastronomy_la_SOURCES = \
    astronomy.F90 \
    include/astronomy_r4.fh \
    include/astronomy_r8.fh \
    include/astronomy.inc

astronomy.$(FC_MODEXT): \
include/astronomy_r4.fh \
include/astronomy_r8.fh \
include/astronomy.inc

BUILT_SOURCES = astronomy_mod.$(FC_MODEXT)
nodist_include_HEADERS = astronomy_mod.$(FC_MODEXT)

EXTRA_DIST = astronomy.tech.ps

include $(top_srcdir)/mkmods.mk
