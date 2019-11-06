# This is an automake file for the monin_obukhov directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libmonin_obukhov.la libmonin_obukhov_kernel.la

libmonin_obukhov_la_SOURCES = monin_obukhov.F90
libmonin_obukhov_kernel_la_SOURCES = monin_obukhov_kernel.F90

# Note that the name of the mod is different from the name of the F90
# code for monin_obukhov_kernel.F90. Also note that the mod file for
# this one does not have "_mod" in the name.
monin_obukhov_mod.mod: monin_obukhov.$(OBJEXT)
monin_obukhov_inter.mod: monin_obukhov_kernel.$(OBJEXT)
monin_obukhov.$(OBJEXT): monin_obukhov_inter.mod

# Mod files are built and then installed as headers.
MODFILES = monin_obukhov_mod.mod monin_obukhov_inter.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
