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

# This is an automake file for the axis_utils directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/axis_utils/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libaxis_utils.la

libaxis_utils_la_SOURCES = \
  axis_utils2.F90 \
  include/axis_utils2_r4.fh \
  include/axis_utils2_r8.fh \
  include/axis_utils2.inc

axis_utils2.$(FC_MODEXT) : \
include/axis_utils2_r4.fh \
include/axis_utils2_r8.fh \
include/axis_utils2.inc

# Mod file depends on its o file, is built and then installed.
nodist_include_HEADERS = axis_utils2_mod.$(FC_MODEXT)
BUILT_SOURCES = axis_utils2_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
