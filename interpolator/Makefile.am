# This is an automake file for the interpolator directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/time_interp
AM_CPPFLAGS += -I${top_builddir}/constants

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libinterpolator.la

# The convenience library depends on its source.
libinterpolator_la_SOURCES = interpolator.F90

# Mod file depends on its o file, is built and then installed.
interpolator_mod.mod: interpolator.$(OBJEXT)
BUILT_SOURCES = interpolator_mod.mod
include_HEADERS = interpolator_mod.mod

CLEANFILES = *.mod
