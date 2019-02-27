# This is an automake file for the time_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/horiz_interp

noinst_LTLIBRARIES = libtime_interp.la libtime_interp_external.la

libtime_interp_la_SOURCES = time_interp.F90
libtime_interp_external_la_SOURCES = time_interp_external.F90

time_interp_mod.mod: time_interp.$(OBJEXT)
time_interp_external_mod.mod: time_interp_external.$(OBJEXT)

MODFILES = time_interp_mod.mod time_interp_external_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
