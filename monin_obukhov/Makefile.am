#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the monin_obukhov directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

noinst_LTLIBRARIES = libmonin_obukhov.la

libmonin_obukhov_la_SOURCES = \
  monin_obukhov.F90 \
  monin_obukhov_inter.F90

monin_obukhov_mod.$(FC_MODEXT): monin_obukhov_inter.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  monin_obukhov_inter.$(FC_MODEXT) \
  monin_obukhov_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = monin_obukhov.tech.ps

include $(top_srcdir)/mkmods.mk
