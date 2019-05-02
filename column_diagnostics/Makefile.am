# This is an automake file for the column_diagnostics directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libcolumn_diagnostics.la

# The convenience library depends on its source.
libcolumn_diagnostics_la_SOURCES = column_diagnostics.F90

# Mod file depends on its o file, is built and then installed.
column_diagnostics_mod.mod: column_diagnostics.$(OBJEXT)
BUILT_SOURCES = column_diagnostics_mod.mod
include_HEADERS = column_diagnostics_mod.mod

CLEANFILES = *.mod
