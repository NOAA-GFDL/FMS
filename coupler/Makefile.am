# This is an automake file for the coupler directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/data_override
AM_CPPFLAGS += -I${top_builddir}/field_manager

noinst_LTLIBRARIES = libcoupler_types.la libensemble_manager.la	\
libatmos_ocean_fluxes.la

libcoupler_types_la_SOURCES = coupler_types.F90
libensemble_manager_la_SOURCES = ensemble_manager.F90
libatmos_ocean_fluxes_la_SOURCES = atmos_ocean_fluxes.F90

CLEANFILES = *.mod
