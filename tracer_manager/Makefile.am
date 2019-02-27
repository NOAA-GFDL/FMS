# This is an automake file for the tracer_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/field_manager

noinst_LTLIBRARIES = libtracer_manager.la

libtracer_manager_la_SOURCES = tracer_manager.F90

tracer_manager_mod.mod: tracer_manager.$(OBJEXT)
BUILT_SOURCES = tracer_manager_mod.mod
include_HEADERS = tracer_manager_mod.mod

CLEANFILES = *.mod
