# This is an automake file for the topography directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/horiz_interp

noinst_LTLIBRARIES = libgaussian_topog.la libtopography.la

libtopography_la_SOURCES = topography.F90
libgaussian_topog_la_SOURCES = gaussian_topog.F90
