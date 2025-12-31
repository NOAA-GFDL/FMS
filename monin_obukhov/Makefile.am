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

# This is an automake file for the monin_obukhov directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/monin_obukhov/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

noinst_LTLIBRARIES = libmonin_obukhov.la

libmonin_obukhov_la_SOURCES = \
        monin_obukhov.F90 \
        monin_obukhov_inter.F90 \
        include/monin_obukhov_r4.fh \
        include/monin_obukhov_r8.fh \
        include/monin_obukhov.inc \
        include/monin_obukhov_inter_r4.fh \
        include/monin_obukhov_inter_r8.fh \
        include/monin_obukhov_inter.inc

monin_obukhov_inter.$(FC_MODEXT): \
        include/monin_obukhov_inter_r4.fh \
        include/monin_obukhov_inter_r8.fh \
        include/monin_obukhov_inter.inc

monin_obukhov_mod.$(FC_MODEXT): \
        monin_obukhov_inter.$(FC_MODEXT) \
        include/monin_obukhov_r4.fh \
        include/monin_obukhov_r8.fh \
        include/monin_obukhov.inc

# Mod files are built and then installed as headers.
MODFILES = \
        monin_obukhov_mod.$(FC_MODEXT) \
        monin_obukhov_inter.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = monin_obukhov.tech.ps

include $(top_srcdir)/mkmods.mk
