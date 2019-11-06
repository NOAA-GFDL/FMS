# This is an automake file for the diag_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdiag_data.la libdiag_axis.la libdiag_grid.la \
libdiag_output.la libdiag_util.la libdiag_manifest.la \
libdiag_table.la libdiag_manager.la

# Each convenience library depends on its source.
libdiag_axis_la_SOURCES = diag_axis.F90
libdiag_data_la_SOURCES = diag_data.F90
libdiag_grid_la_SOURCES = diag_grid.F90
libdiag_manager_la_SOURCES = diag_manager.F90
libdiag_manifest_la_SOURCES = diag_manifest.F90
libdiag_output_la_SOURCES = diag_output.F90
libdiag_table_la_SOURCES = diag_table.F90
libdiag_util_la_SOURCES = diag_util.F90

# Each mod file depends on the .o file.
diag_axis_mod.mod: diag_axis.$(OBJEXT)
diag_grid_mod.mod: diag_grid.$(OBJEXT)
diag_manifest_mod.mod: diag_manifest.$(OBJEXT)
diag_table_mod.mod: diag_table.$(OBJEXT)
diag_data_mod.mod: diag_data.$(OBJEXT)
diag_manager_mod.mod: diag_manager.$(OBJEXT)
diag_output_mod.mod: diag_output.$(OBJEXT)
diag_util_mod.mod: diag_util.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
diag_axis.$(OBJEXT): diag_data_mod.mod
diag_manifest.$(OBJEXT): diag_data_mod.mod
diag_output.$(OBJEXT): diag_axis_mod.mod diag_data_mod.mod
diag_util.$(OBJEXT): diag_data_mod.mod diag_axis_mod.mod diag_output_mod.mod diag_grid_mod.mod
diag_table.$(OBJEXT): diag_data_mod.mod diag_util_mod.mod
diag_manager.$(OBJEXT): diag_axis_mod.mod diag_data_mod.mod diag_util_mod.mod diag_output_mod.mod diag_grid_mod.mod diag_table_mod.mod diag_manifest_mod.mod

# Mod files are built and then installed as headers.
MODFILES = diag_data_mod.mod diag_axis_mod.mod diag_grid_mod.mod \
diag_output_mod.mod diag_util_mod.mod diag_manifest_mod.mod \
diag_table_mod.mod diag_manager_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
