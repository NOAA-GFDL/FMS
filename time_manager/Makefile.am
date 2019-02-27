# This is an automake file for the time_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp

noinst_LTLIBRARIES = libtime_manager.la libget_cal_time.la

libtime_manager_la_SOURCES = time_manager.F90
libget_cal_time_la_SOURCES = get_cal_time.F90

include_HEADERS = get_cal_time_mod.mod time_manager_mod.mod

CLEANFILES = *.mod
