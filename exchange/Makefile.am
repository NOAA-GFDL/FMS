# This is an automake file for the exchange directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mosaic

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libstock_constants.la libxgrid.la

# Each convenience library depends on its source.
libstock_constants_la_SOURCES = stock_constants.F90
libxgrid_la_SOURCES = xgrid.F90 stock_constants_mod.mod

# Each mod file depends on the .o file.
stock_constants_mod.mod: stock_constants.$(OBJEXT)
xgrid_mod.mod: xgrid.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
xgrid.$(OBJEXT): stock_constants_mod.mod

# Mod files are built and then installed as headers.
MODFILES = stock_constants_mod.mod xgrid_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
