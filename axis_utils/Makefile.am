# This is an automake file for the axis_utils directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms2_io

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libaxis_utils.la libaxis_utils2.la

# The convenience library depends on its source.
libaxis_utils_la_SOURCES = axis_utils.F90
libaxis_utils2_la_SOURCES = axis_utils2.F90

# Mod file depends on its o file, is built and then installed.
axis_utils_mod.mod: axis_utils.$(OBJEXT)
axis_utils2_mod.mod: axis_utils2.$(OBJEXT)
include_HEADERS = axis_utils_mod.mod axis_utils2_mod.mod
BUILT_SOURCES = axis_utils_mod.mod axis_utils2_mod.mod

CLEANFILES = *.mod
