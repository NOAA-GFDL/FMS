# This is an automake file for the memutils directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

# The convenience library depends on its source.
libmemutils_la_SOURCES = memutils.F90
if APPLE
libmemuse_c_la_SOURCES = memuse.c
noinst_LTLIBRARIES = libmemutils.la libmemuse_c.la
else
noinst_LTLIBRARIES = libmemutils.la
endif

# Mod file depends on its o file, is built and then installed.
memutils_mod.mod: memutils.$(OBJEXT)
BUILT_SOURCES = memutils_mod.mod
include_HEADERS = memutils_mod.mod

CLEANFILES = *.mod
