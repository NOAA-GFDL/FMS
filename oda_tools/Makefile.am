# This is an automake file for the oda_tools directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/field_manager

noinst_LTLIBRARIES = liboda_types.la libwrite_ocean_data.la	\
liboda_core.la

liboda_types_la_SOURCES = oda_types.F90
#liboda_core_ecda_la_SOURCES = oda_core_ecda.F90
liboda_core_la_SOURCES = oda_core.F90
libwrite_ocean_data_la_SOURCES = write_ocean_data.F90

oda_types_mod.mod: oda_types.$(OBJEXT)
write_ocean_data_mod.mod: write_ocean_data.$(OBJEXT)
oda_core_mod.mod: oda_core.$(OBJEXT)

MODFILES = oda_types_mod.mod write_ocean_data_mod.mod	\
oda_core_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
