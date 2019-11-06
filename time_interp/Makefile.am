# This is an automake file for the time_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/horiz_interp

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtime_interp.la libtime_interp_external.la

# Each convenience library depends on its source.
libtime_interp_la_SOURCES = time_interp.F90
libtime_interp_external_la_SOURCES = time_interp_external.F90

# Each mod file depends on the .o file.
time_interp_mod.mod: time_interp.$(OBJEXT)
time_interp_external_mod.mod: time_interp_external.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
time_interp_external.$(OBJEXT):  time_interp_mod.mod

# Mod files are built and then installed as headers.
MODFILES = time_interp_mod.mod time_interp_external_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
