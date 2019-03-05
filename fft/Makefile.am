# This is an automake file for the fft directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfft99.la libfft.la

# Each convenience library depends on its source.
libfft_la_SOURCES = fft.F90
libfft99_la_SOURCES = fft99.F90

# Each mod file depends on the .o file.
fft99_mod.mod: fft99.$(OBJEXT)
fft_mod.mod: fft.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
fft.$(OBJEXT): fft99_mod.mod

# Mod files are built and then installed as headers.
MODFILES = fft99_mod.mod fft_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
