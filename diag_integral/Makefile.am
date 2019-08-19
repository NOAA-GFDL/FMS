# This is an automake file for the diag_integral directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/coupler

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libdiag_integral.la

# The convenience library depends on its source.
libdiag_integral_la_SOURCES = diag_integral.F90

# Mod file depends on its o file, is built and then installed.
diag_integral_mod.mod: diag_integral.$(OBJEXT)
include_HEADERS = diag_integral_mod.mod
BUILT_SOURCES = diag_integral_mod.mod

EXTRA_DIST = diag_integral.html

CLEANFILES = *.mod
