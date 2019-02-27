# This is an automake file for the memutils directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

noinst_LTLIBRARIES = libmemutils.la

libmemutils_la_SOURCES = memutils.F90

memutils_mod.mod: memutils.$(OBJEXT)
BUILT_SOURCES = memutils_mod.mod
include_HEADERS = memutils_mod.mod

CLEANFILES = *.mod
