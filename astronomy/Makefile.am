# This is an automake file for the astronomy directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libastronomy.la

# The convenience library depends on its source.
libastronomy_la_SOURCES = astronomy.F90

# Mod file depends on its o file, is built and then installed.
astronomy_mod.mod: astronomy.$(OBJEXT)
BUILT_SOURCES = astronomy_mod.mod
include_HEADERS = astronomy_mod.mod

EXTRA_DIST = astronomy.tech.ps

CLEANFILES = *.mod
