# This is an automake file for the fft directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp

noinst_LTLIBRARIES = libfft99.la libfft.la

libfft_la_SOURCES = fft.F90
libfft99_la_SOURCES = fft99.F90

CLEANFILES = *.mod
