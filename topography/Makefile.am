# This is an automake file for the topography directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp

noinst_LTLIBRARIES = libgaussian_topog.la libtopography.la

libtopography_la_SOURCES = topography.F90
libgaussian_topog_la_SOURCES = gaussian_topog.F90

gaussian_topog_mod.mod: gaussian_topog.$(OBJEXT)
topography_mod.mod: topography.$(OBJEXT)

MODFILES = gaussian_topog_mod.mod topography_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
