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
  drifters_compute_k.fh \
  drifters_push.fh \
  drifters_set_field.fh \
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
