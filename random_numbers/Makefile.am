# This is an automake file for the random_numbers directory of the FMS
# package.

# Ed Hartnett 2/28/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/time_manager

# Build these uninstalled convenience library.
noinst_LTLIBRARIES = libMersenneTwister.la librandom_numbers.la

# Each convenience library depends on its source.
librandom_numbers_la_SOURCES = random_numbers.F90
libMersenneTwister_la_SOURCES = MersenneTwister.F90

# Each mod file depends on the .o file.
random_numbers_mod.mod: random_numbers.$(OBJEXT)
mersennetwister_mod.mod: MersenneTwister.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
random_numbers.$(OBJEXT): mersennetwister_mod.mod

# Mod files are built and then installed as headers.
MODFILES = mersennetwister_mod.mod random_numbers_mod.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

CLEANFILES = *.mod
