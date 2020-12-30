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
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)


# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdrifters.la

# Each convenience library depends on its source.
libdrifters_la_SOURCES = drifters.F90 \
  drifters_core.F90 \
  drifters_comm.F90 \
  drifters_input.F90 \
  drifters_io.F90 \
  cloud_interpolator.F90 \
  quicksort.F90 \
  drifters_compute_k.h \
  drifters_push.h \
  drifters_set_field.h \
  fms_switches.h

drifters_mod.$(FC_MODEXT): drifters_core_mod.$(FC_MODEXT) drifters_input_mod.$(FC_MODEXT) drifters_io_mod.$(FC_MODEXT) \
                  drifters_comm_mod.$(FC_MODEXT) cloud_interpolator_mod.$(FC_MODEXT)
drifters_comm_mod.$(FC_MODEXT): drifters_core_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  drifters_core_mod.$(FC_MODEXT) \
  drifters_input_mod.$(FC_MODEXT) \
  drifters_comm_mod.$(FC_MODEXT) \
  drifters_io_mod.$(FC_MODEXT) \
  cloud_interpolator_mod.$(FC_MODEXT) \
  drifters_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES) $(DRIFTERS_HDRS)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = README.md

include $(top_srcdir)/mkmods.mk
