# This is an automake file for the fms directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/memutils

noinst_LTLIBRARIES = libfms_io.la libfms.la

fms_io_mod.mod: fms_io.$(OBJEXT)
fms_mod.mod: fms.$(OBJEXT)
BUILT_SOURCES = fms_io_mod.mod fms_mod.mod
include_HEADERS = fms_io_mod.mod fms_mod.mod

libfms_la_SOURCES = fms.F90
libfms_io_la_SOURCES = fms_io.F90

CLEANFILES = *.mod
