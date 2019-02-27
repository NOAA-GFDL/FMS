# This is an automake file for the platform directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include

noinst_LTLIBRARIES = libplatform.la

libplatform_la_SOURCES = platform.F90

platform_mod.mod: platform.$(OBJEXT)
BUILT_SOURCES = platform_mod.mod
include_HEADERS = platform_mod.mod

CLEANFILES = *.mod
