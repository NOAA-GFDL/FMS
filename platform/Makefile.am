# This is an automake file for the platform directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h files.
AM_CPPFLAGS = -I${top_srcdir}/include

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libplatform.la

# The convenience library depends on its source.
libplatform_la_SOURCES = platform.F90

# Mod file depends on its o file, is built and then installed.
platform_mod.mod: platform.$(OBJEXT)
BUILT_SOURCES = platform_mod.mod
include_HEADERS = platform_mod.mod

CLEANFILES = *.mod
