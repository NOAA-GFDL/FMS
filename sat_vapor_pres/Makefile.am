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

# This is an automake file for the sat_vapor_pres directory of the FMS
# package.

# Ed Hartnett 2/28/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/sat_vapor_pres/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience library.
noinst_LTLIBRARIES = libsat_vapor_pres.la

# Each convenience library depends on its source.
libsat_vapor_pres_la_SOURCES = \
	sat_vapor_pres.F90 \
	include/sat_vapor_pres_r4.fh \
	include/sat_vapor_pres_r8.fh \
	include/sat_vapor_pres.inc \
	sat_vapor_pres_k.F90 \
	include/sat_vapor_pres_k_r4.fh \
	include/sat_vapor_pres_k_r8.fh \
	include/sat_vapor_pres_k.inc

# Some mods are dependant on other mods in this dir.
sat_vapor_pres_mod.$(FC_MODEXT): \
	sat_vapor_pres_k_mod.$(FC_MODEXT) \
	include/sat_vapor_pres_r4.fh \
	include/sat_vapor_pres_r8.fh \
	include/sat_vapor_pres.inc
sat_vapor_pres_k_mod.$(FC_MODEXT): \
	include/sat_vapor_pres_k_r4.fh \
	include/sat_vapor_pres_k_r8.fh \
	include/sat_vapor_pres_k.inc

# Mod files are built and then installed as headers.
MODFILES = \
	sat_vapor_pres_k_mod.$(FC_MODEXT) \
	sat_vapor_pres_mod.$(FC_MODEXT)
BUILT_SOURCES = $(MODFILES)
nodist_include_HEADERS = $(MODFILES)

include $(top_srcdir)/mkmods.mk
