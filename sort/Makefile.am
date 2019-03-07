# This is an automake file for the sort directory of the FMS
# package.

# Ed Hartnett 3/7/19

# Include .h files.
AM_CPPFLAGS = -I${top_srcdir}/include

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libsort_quicksort.la

# The convenience library depends on its source.
libsort_quicksort_la_SOURCES = sort_quicksort.f90

# Mod file depends on its o file, is built and then installed.
quicksort_mod.mod: sort_quicksort.$(OBJEXT)
BUILT_SOURCES = quicksort_mod.mod
include_HEADERS = quicksort_mod.mod

CLEANFILES = *.mod
