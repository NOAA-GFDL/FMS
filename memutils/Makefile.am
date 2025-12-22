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

# This is an automake file for the memutils directory of the FMS package.

# Copyright 2019-2020, GFDL, Ed Hartnett, Seth Underwood

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# The convenience library depends on its source.
noinst_LTLIBRARIES = libmemutils.la
libmemutils_la_SOURCES = memutils.F90

BUILT_SOURCES = memutils_mod.$(FC_MODEXT)
nodist_include_HEADERS = memutils_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
