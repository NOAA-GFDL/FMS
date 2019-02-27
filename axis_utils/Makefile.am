# This is an automake file for the axis_utils directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libaxis_utils.la

# The convenience library depends on its source.
libaxis_utils_la_SOURCES = axis_utils.F90

# Mod file depends on its o file, is built and then installed.
axis_utils_mod.mod: axis_utils.$(OBJEXT)
include_HEADERS = axis_utils_mod.mod
BUILT_SOURCES = axis_utils_mod.mod

CLEANFILES = *.mod
