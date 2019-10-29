# This is an automake file for the memutils directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libfms_affinity.la libfms_affinity_c.la

# The convenience library depends on its source.
libfms_affinity_la_SOURCES = fms_affinity.F90
libfms_affinity_c_la_SOURCES = affinity.c

# Mod file depends on its .o file, is built and then installed.
fms_affinity_mod.mod: fms_affinity.$(OBJEXT)

# Mod files are built and then installed as headers
MODFILES = fms_affinity_mod.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

CLEANFILES = *.mod
