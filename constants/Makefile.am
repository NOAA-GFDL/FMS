# This is an automake file for the constants directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/platform

noinst_LTLIBRARIES = libconstants.la

libconstants_la_SOURCES = constants.F90

constants_mod.mod: constants.$(OBJEXT)
BUILT_SOURCES = constants_mod.mod
include_HEADERS = constants_mod.mod

CLEANFILES = *.mod
