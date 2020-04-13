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

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libmonin_obukhov.la libmonin_obukhov_kernel.la

libmonin_obukhov_la_SOURCES = monin_obukhov.F90
libmonin_obukhov_kernel_la_SOURCES = monin_obukhov_kernel.F90

# Note that the name of the mod is different from the name of the F90
# code for monin_obukhov_kernel.F90. Also note that the mod file for
# this one does not have "_mod" in the name.
monin_obukhov_mod.mod: monin_obukhov.lo
monin_obukhov_inter.mod: monin_obukhov_kernel.lo
monin_obukhov.lo: monin_obukhov_inter.mod

# Mod files are built and then installed as headers.
MODFILES = monin_obukhov_mod.mod monin_obukhov_inter.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
