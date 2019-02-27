# This is an automake file for the tracer_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/field_manager

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libtracer_manager.la

# The convenience library depends on its source.
libtracer_manager_la_SOURCES = tracer_manager.F90

# Mod file depends on its o file, is built and then installed.
tracer_manager_mod.mod: tracer_manager.$(OBJEXT)
BUILT_SOURCES = tracer_manager_mod.mod
include_HEADERS = tracer_manager_mod.mod

CLEANFILES = *.mod
