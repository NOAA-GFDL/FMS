# This is an automake file for the block_control directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libblock_control.la

# The convenience library depends on its source.
libblock_control_la_SOURCES = block_control.F90

# Mod file depends on its o file, is built and then installed.
block_control_mod.mod: block_control.$(OBJEXT)
BUILT_SOURCES = block_control_mod.mod
include_HEADERS = block_control_mod.mod

CLEANFILES = *.mod
