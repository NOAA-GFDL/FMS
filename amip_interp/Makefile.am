# This is an automake file for the amip_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_interp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/platform

noinst_LTLIBRARIES = libamip_interp.la

libamip_interp_la_SOURCES = amip_interp.F90

amip_interp_mod.mod: amip_interp.$(OBJEXT)
BUILT_SOURCES = amip_interp_mod.mod
include_HEADERS = amip_interp_mod.mod

CLEANFILES = *.mod
