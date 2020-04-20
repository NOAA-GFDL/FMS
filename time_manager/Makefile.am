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

# This is an automake file for the time_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtime_manager.la libget_cal_time.la

# Each convenience library depends on its source.
libtime_manager_la_SOURCES = time_manager.F90
libget_cal_time_la_SOURCES = get_cal_time.F90

# Each mod file depends on the .o file.
time_manager_mod.mod: time_manager.lo
get_cal_time_mod.mod: get_cal_time.lo

# Some mods are dependant on other mods in this dir.
get_cal_time.lo: time_manager_mod.mod

# Mod files are built and then installed as headers.
MODFILES = time_manager_mod.mod get_cal_time_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
