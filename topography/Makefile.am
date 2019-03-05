# This is an automake file for the topography directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libgaussian_topog.la libtopography.la

# Each convenience library depends on its source.
libtopography_la_SOURCES = topography.F90
libgaussian_topog_la_SOURCES = gaussian_topog.F90

# Each mod file depends on the .o file.
gaussian_topog_mod.mod: gaussian_topog.$(OBJEXT)
topography_mod.mod: topography.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
topography.$(OBJEXT): gaussian_topog_mod.mod

# Mod files are built and then installed as headers.
MODFILES = gaussian_topog_mod.mod topography_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
