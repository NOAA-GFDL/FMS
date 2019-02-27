# This is an automake file for the station_data directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/time_manager

noinst_LTLIBRARIES = libstation_data.la

libstation_data_la_SOURCES = station_data.F90

station_data_mod.mod: station_data.$(OBJEXT)
BUILT_SOURCES = station_data_mod.mod
include_HEADERS = station_data_mod.mod

CLEANFILES = *.mod
