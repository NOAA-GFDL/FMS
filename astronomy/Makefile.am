# This is an automake file for the astronomy directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp

noinst_LTLIBRARIES = libastronomy.la

libastronomy_la_SOURCES = astronomy.F90

EXTRA_DIST = astronomy.tech.ps
