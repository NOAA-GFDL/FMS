# This is an automake file for the constants directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/platform

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libconstants.la

# The convenience library depends on its source.
libconstants_la_SOURCES = constants.F90

# Mod file depends on its o file, is built and then installed.
constants_mod.mod: constants.$(OBJEXT)
BUILT_SOURCES = constants_mod.mod
include_HEADERS = constants_mod.mod

CLEANFILES = *.mod
