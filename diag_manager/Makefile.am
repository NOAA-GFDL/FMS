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
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdiag_manager.la

# Each convenience library depends on its source.
libdiag_manager_la_SOURCES = \
  diag_axis.F90 \
  diag_data.F90 \
  diag_grid.F90 \
  diag_manager.F90 \
  diag_output.F90 \
  diag_table.F90 \
  diag_util.F90

# Some mods are dependant on other mods in this dir.
diag_axis_mod.$(FC_MODEXT): diag_data_mod.$(FC_MODEXT)
diag_output_mod.$(FC_MODEXT): diag_axis_mod.$(FC_MODEXT) diag_data_mod.$(FC_MODEXT)
diag_util_mod.$(FC_MODEXT): diag_data_mod.$(FC_MODEXT) diag_axis_mod.$(FC_MODEXT) diag_output_mod.$(FC_MODEXT) diag_grid_mod.$(FC_MODEXT)
diag_table_mod.$(FC_MODEXT): diag_data_mod.$(FC_MODEXT) diag_util_mod.$(FC_MODEXT)
diag_manager_mod.$(FC_MODEXT): diag_axis_mod.$(FC_MODEXT) diag_data_mod.$(FC_MODEXT) diag_util_mod.$(FC_MODEXT) \
	                   diag_output_mod.$(FC_MODEXT) diag_grid_mod.$(FC_MODEXT) diag_table_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  diag_data_mod.$(FC_MODEXT) \
  diag_axis_mod.$(FC_MODEXT) \
  diag_grid_mod.$(FC_MODEXT) \
  diag_output_mod.$(FC_MODEXT) \
  diag_util_mod.$(FC_MODEXT) \
  diag_table_mod.$(FC_MODEXT) \
  diag_manager_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

include $(top_srcdir)/mkmods.mk
