# This is an automake file for the coupler directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/data_override
AM_CPPFLAGS += -I${top_builddir}/field_manager

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcoupler_types.la libensemble_manager.la \
libatmos_ocean_fluxes.la

# Each convenience library depends on its source.
libcoupler_types_la_SOURCES = coupler_types.F90
libensemble_manager_la_SOURCES = ensemble_manager.F90
libatmos_ocean_fluxes_la_SOURCES = atmos_ocean_fluxes.F90

# Each mod file depends on the .o file.
coupler_types_mod.mod: coupler_types.$(OBJEXT)
ensemble_manager_mod.mod: ensemble_manager.$(OBJEXT)
atmos_ocean_fluxes_mod.mod: atmos_ocean_fluxes.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
atmos_ocean_fluxes.$(OBJEXT): coupler_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = coupler_types_mod.mod ensemble_manager_mod.mod \
atmos_ocean_fluxes_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
