# This is an automake file for the horiz_interp directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libhoriz_interp_type.la \
libhoriz_interp_bicubic.la libhoriz_interp_bilinear.la \
libhoriz_interp_conserve.la libhoriz_interp_spherical.la \
libhoriz_interp.la

# Each convenience library depends on its source.
libhoriz_interp_bicubic_la_SOURCES = horiz_interp_bicubic.F90
libhoriz_interp_bilinear_la_SOURCES = horiz_interp_bilinear.F90
libhoriz_interp_conserve_la_SOURCES = horiz_interp_conserve.F90
libhoriz_interp_la_SOURCES = horiz_interp.F90
libhoriz_interp_spherical_la_SOURCES = horiz_interp_spherical.F90
libhoriz_interp_type_la_SOURCES = horiz_interp_type.F90

# Each mod file depends on the .o file.
horiz_interp_type_mod.mod: horiz_interp_type.$(OBJEXT)
horiz_interp_bicubic_mod.mod: horiz_interp_bicubic.$(OBJEXT)
horiz_interp_bilinear_mod.mod: horiz_interp_bilinear.$(OBJEXT)
horiz_interp_conserve_mod.mod: horiz_interp_conserve.$(OBJEXT)
horiz_interp_spherical_mod.mod: horiz_interp_spherical.$(OBJEXT)
horiz_interp_mod.mod: horiz_interp.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
horiz_interp_bicubic.$(OBJEXT): horiz_interp_type_mod.mod
horiz_interp_bilinear.$(OBJEXT): horiz_interp_type_mod.mod
horiz_interp_conserve.$(OBJEXT): horiz_interp_type_mod.mod
horiz_interp_spherical.$(OBJEXT): horiz_interp_type_mod.mod
horiz_interp.$(OBJEXT): horiz_interp_bicubic_mod.mod horiz_interp_type_mod.mod \
  horiz_interp_bilinear_mod.mod horiz_interp_conserve_mod.mod horiz_interp_spherical_mod.mod

# Mod files are built and then installed as headers.
MODFILES = horiz_interp_type_mod.mod horiz_interp_bicubic_mod.mod \
horiz_interp_bilinear_mod.mod horiz_interp_conserve_mod.mod \
horiz_interp_spherical_mod.mod horiz_interp_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
