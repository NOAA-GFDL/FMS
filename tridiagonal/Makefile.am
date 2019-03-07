# This is an automake file for the tridiagonal directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/platform

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libtridiagonal.la

# The convenience library depends on its source.
libtridiagonal_la_SOURCES = tridiagonal.F90

# Mod file depends on its o file, is built and then installed.
tridiagonal_mod.mod: tridiagonal.$(OBJEXT)
BUILT_SOURCES = tridiagonal_mod.mod
include_HEADERS = tridiagonal_mod.mod

CLEANFILES = *.mod
