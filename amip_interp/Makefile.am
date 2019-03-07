# This is an automake file for the amip_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_interp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/platform

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libamip_interp.la

# The convenience library depends on its source.
libamip_interp_la_SOURCES = amip_interp.F90

# Mod file depends on its o file, is built and then installed.
amip_interp_mod.mod: amip_interp.$(OBJEXT)
BUILT_SOURCES = amip_interp_mod.mod
include_HEADERS = amip_interp_mod.mod

CLEANFILES = *.mod
