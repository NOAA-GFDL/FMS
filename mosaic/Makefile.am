# This is an automake file for the mosaic directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms

noinst_LTLIBRARIES = libmosaic.la libgrid.la libgradient.la

libmosaic_la_SOURCES = mosaic.F90
libgrid_la_SOURCES = grid.F90
libgradient_la_SOURCES = gradient.F90
