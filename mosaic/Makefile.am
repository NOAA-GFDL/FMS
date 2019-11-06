# This is an automake file for the mosaic directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmosaic.la libgrid.la libgradient.la \
libmosaic_c.la

# Each convenience library depends on its source. Lump all the C code
# into libmosaic_c.la.
libmosaic_la_SOURCES = mosaic.F90
libgrid_la_SOURCES = grid.F90
libgradient_la_SOURCES = gradient.F90
libmosaic_c_la_SOURCES = create_xgrid.c gradient_c2l.c interp.c \
mosaic_util.c read_mosaic.c constant.h create_xgrid.h gradient_c2l.h \
interp.h mosaic_util.h read_mosaic.h

# Each mod file depends on the .o file.
mosaic_mod.mod: mosaic.$(OBJEXT)
grid_mod.mod: grid.$(OBJEXT)
gradient_mod.mod: gradient.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
grid.$(OBJEXT): mosaic_mod.mod

# Mod files are built and then installed as headers.
MODFILES = mosaic_mod.mod grid_mod.mod gradient_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
