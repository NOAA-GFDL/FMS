# This is an automake file for the column_diagnostics directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libcolumn_diagnostics.la

libcolumn_diagnostics_la_SOURCES = column_diagnostics.F90

column_diagnostics_mod.mod: column_diagnostics.$(OBJEXT)
BUILT_SOURCES = column_diagnostics_mod.mod
include_HEADERS = column_diagnostics_mod.mod

CLEANFILES = *.mod
