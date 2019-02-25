# This is an automake file for the field_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms

noinst_LTLIBRARIES = libfield_manager.la

libfield_manager_la_SOURCES = field_manager.F90 fm_util.F90

EXTRA_DIST = parse.inc
