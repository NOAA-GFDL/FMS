# This is an automake file for the time_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtime_manager.la libget_cal_time.la

# Each convenience library depends on its source.
libtime_manager_la_SOURCES = time_manager.F90
libget_cal_time_la_SOURCES = get_cal_time.F90

# Each mod file depends on the .o file.
time_manager_mod.mod: time_manager.$(OBJEXT)
get_cal_time_mod.mod: get_cal_time.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
get_cal_time.$(OBJEXT): time_manager_mod.mod

# Mod files are built and then installed as headers.
MODFILES = time_manager_mod.mod get_cal_time_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
