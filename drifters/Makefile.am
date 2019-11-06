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
drifters_core_mod.mod: drifters_core.$(OBJEXT)
drifters_input_mod.mod: drifters_input.$(OBJEXT)
drifters_comm_mod.mod: drifters_comm.$(OBJEXT)
drifters_io_mod.mod: drifters_io.$(OBJEXT)
cloud_interpolator_mod.mod: cloud_interpolator.$(OBJEXT)
drifters_mod.mod: drifters.$(OBJEXT)

# Mod files are built and then installed as headers.
MODFILES = drifters_core_mod.mod drifters_input_mod.mod \
drifters_comm_mod.mod drifters_io_mod.mod cloud_interpolator_mod.mod \
drifters_mod.mod
include_HEADERS = $(MODFILES) $(DRIFTERS_HDRS)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod *__genmod.f90
