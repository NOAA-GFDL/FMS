# This is an automake file for the diag_integral directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libdiag_integral.la

libdiag_integral_la_SOURCES = diag_integral.F90

EXTRA_DIST = diag_integral.html
