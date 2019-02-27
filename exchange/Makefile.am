# This is an automake file for the exchange directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mosaic

noinst_LTLIBRARIES = libstock_constants.la libxgrid.la

libstock_constants_la_SOURCES = stock_constants.F90
libxgrid_la_SOURCES = xgrid.F90

stock_constants_mod.mod: stock_constants.$(OBJEXT)
xgrid_mod.mod: xgrid.$(OBJEXT)

MODFILES = stock_constants_mod.mod xgrid_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
