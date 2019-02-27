# This is an automake file for the monin_obukhov directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms

# Allow constants.mod to be found.
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libmonin_obukhov.la libmonin_obukhov_kernel.la

libmonin_obukhov_la_SOURCES = monin_obukhov.F90
libmonin_obukhov_kernel_la_SOURCES = monin_obukhov_kernel.F90

CLEANFILES = *.mod
