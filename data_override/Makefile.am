# This is an automake file for the data_override directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/time_interp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/axis_utils

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libdata_override.la

# The convenience library depends on its source.
libdata_override_la_SOURCES = data_override.F90

# Mod file depends on its o file, is built and then installed.
data_override_mod.mod: data_override.$(OBJEXT)
BUILT_SOURCES = data_override_mod.mod
include_HEADERS = data_override_mod.mod

CLEANFILES = *.mod
