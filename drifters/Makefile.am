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

# This is an automake file for the drifters directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdrifters_core.la libdrifters_input.la \
libdrifters_comm.la libdrifters_io.la libcloud_interpolator.la \
libdrifters.la libquicksort.la

# Each convenience library depends on its source.
libdrifters_la_SOURCES = drifters.F90
libdrifters_core_la_SOURCES = drifters_core.F90
libdrifters_comm_la_SOURCES = drifters_comm.F90
libdrifters_input_la_SOURCES = drifters_input.F90
libdrifters_io_la_SOURCES = drifters_io.F90
libcloud_interpolator_la_SOURCES = cloud_interpolator.F90
libquicksort_la_SOURCES = quicksort.F90

# These headers are in the drifters directory.
DRIFTERS_HDRS = drifters_compute_k.h drifters_push.h \
drifters_set_field.h fms_switches.h

# Each mod file depends on the .o file.
drifters_core_mod.mod: drifters_core.lo
drifters_input_mod.mod: drifters_input.lo
drifters_comm_mod.mod: drifters_comm.lo
drifters_io_mod.mod: drifters_io.lo
cloud_interpolator_mod.mod: cloud_interpolator.lo
drifters_mod.mod: drifters.lo

# Mod files are built and then installed as headers.
MODFILES = drifters_core_mod.mod drifters_input_mod.mod \
drifters_comm_mod.mod drifters_io_mod.mod cloud_interpolator_mod.mod \
drifters_mod.mod
include_HEADERS = $(MODFILES) $(DRIFTERS_HDRS)
BUILT_SOURCES = $(MODFILES)
CLEANFILES = *.mod *__genmod.f90

