# This is an automake file for the field_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfield_manager.la libfm_util.la

# Each convenience library depends on its source.
libfield_manager_la_SOURCES = field_manager.F90
libfm_util_la_SOURCES = fm_util.F90

# Each mod file depends on the .o file.
field_manager_mod.mod: field_manager.$(OBJEXT)
fm_util_mod.mod: fm_util.$(OBJEXT)

# Mod files are built and then installed as headers.
MODFILES = field_manager_mod.mod fm_util_mod.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

EXTRA_DIST = parse.inc

CLEANFILES = *.mod
