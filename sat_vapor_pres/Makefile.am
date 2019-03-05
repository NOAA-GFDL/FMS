# This is an automake file for the sat_vapor_pres directory of the FMS
# package.

# Ed Hartnett 2/28/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build these uninstalled convenience library.
noinst_LTLIBRARIES = libsat_vapor_pres_k.la libsat_vapor_pres.la

# Each convenience library depends on its source.
libsat_vapor_pres_la_SOURCES = sat_vapor_pres.F90
libsat_vapor_pres_k_la_SOURCES = sat_vapor_pres_k.F90

# Each mod file depends on the .o file.
sat_vapor_pres_mod.mod: sat_vapor_pres.$(OBJEXT)
sat_vapor_pres_k_mod.mod: sat_vapor_pres_k.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
sat_vapor_pres.$(OBJEXT): sat_vapor_pres_k_mod.mod

# Mod files are built and then installed as headers.
MODFILES = sat_vapor_pres_k_mod.mod sat_vapor_pres_mod.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

CLEANFILES = *.mod
