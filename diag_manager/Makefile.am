#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the diag_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/fms2_io
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/platform


# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdiag_data.la libdiag_axis.la libdiag_grid.la	\
libdiag_output.la libdiag_util.la libdiag_table.la libdiag_manager.la

# Each convenience library depends on its source.
libdiag_axis_la_SOURCES = diag_axis.F90
libdiag_data_la_SOURCES = diag_data.F90
libdiag_grid_la_SOURCES = diag_grid.F90
libdiag_manager_la_SOURCES = diag_manager.F90
libdiag_output_la_SOURCES = diag_output.F90
libdiag_table_la_SOURCES = diag_table.F90
libdiag_util_la_SOURCES = diag_util.F90

# Each mod file depends on the .o file.
diag_axis_mod.mod: diag_axis.lo
diag_grid_mod.mod: diag_grid.lo
diag_table_mod.mod: diag_table.lo
diag_data_mod.mod: diag_data.lo
diag_manager_mod.mod: diag_manager.lo
diag_output_mod.mod: diag_output.lo
diag_util_mod.mod: diag_util.lo

# Some mods are dependant on other mods in this dir.
diag_axis.lo: diag_data_mod.mod
diag_output.lo: diag_axis_mod.mod diag_data_mod.mod
diag_util.lo: diag_data_mod.mod diag_axis_mod.mod diag_output_mod.mod diag_grid_mod.mod
diag_table.lo: diag_data_mod.mod diag_util_mod.mod
diag_manager.lo: diag_axis_mod.mod diag_data_mod.mod diag_util_mod.mod diag_output_mod.mod diag_grid_mod.mod diag_table_mod.mod

# Mod files are built and then installed as headers.
MODFILES = diag_data_mod.mod diag_axis_mod.mod diag_grid_mod.mod	\
diag_output_mod.mod diag_util_mod.mod diag_table_mod.mod		\
diag_manager_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
