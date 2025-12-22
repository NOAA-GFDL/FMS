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

# This is an automake file for the random_numbers directory of the FMS
# package.

# Ed Hartnett 2/28/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/random_numbers/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience library.
noinst_LTLIBRARIES = librandom_numbers.la

# Each convenience library depends on its source.
librandom_numbers_la_SOURCES = \
  mersennetwister.F90 \
  random_numbers.F90 \
  include/random_numbers.inc \
  include/random_numbers_r4.fh \
  include/random_numbers_r8.fh

# Some mods are dependant on other mods in this dir.
random_numbers_mod.$(FC_MODEXT): mersennetwister_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  mersennetwister_mod.$(FC_MODEXT) \
  random_numbers_mod.$(FC_MODEXT)
BUILT_SOURCES = $(MODFILES)
nodist_include_HEADERS = $(MODFILES)

include $(top_srcdir)/mkmods.mk
