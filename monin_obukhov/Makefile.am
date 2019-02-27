# This is an automake file for the monin_obukhov directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms

# Allow constants.mod to be found.
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libmonin_obukhov.la libmonin_obukhov_kernel.la

libmonin_obukhov_la_SOURCES = monin_obukhov.F90
libmonin_obukhov_kernel_la_SOURCES = monin_obukhov_kernel.F90

monin_obukhov_mod.mod: monin_obukhov.$(OBJEXT)
monin_obukhov_kernel_mod.mod: monin_obukhov_kernel.$(OBJEXT)

MODFILES = monin_obukhov_mod.mod monin_obukhov_kernel_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
