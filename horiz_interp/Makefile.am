# This is an automake file for the horiz_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libhoriz_interp_type.la			\
libhoriz_interp_bicubic.la libhoriz_interp_bilinear.la		\
libhoriz_interp_conserve.la libhoriz_interp_spherical.la	\
libhoriz_interp.la

libhoriz_interp_bicubic_la_SOURCES = horiz_interp_bicubic.F90
libhoriz_interp_bilinear_la_SOURCES = horiz_interp_bilinear.F90
libhoriz_interp_conserve_la_SOURCES = horiz_interp_conserve.F90
libhoriz_interp_la_SOURCES = horiz_interp.F90
libhoriz_interp_spherical_la_SOURCES = horiz_interp_spherical.F90
libhoriz_interp_type_la_SOURCES = horiz_interp_type.F90

CLEANFILES = *.mod
