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

# This is an automake file for the field_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/field_manager/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfield_manager.la

# Each convenience library depends on its source.
libfield_manager_la_SOURCES = \
	field_manager.F90 \
	fm_yaml.F90 \
	fm_util.F90 \
	parse.inc \
	include/field_manager.inc \
	include/field_manager_r4.fh \
	include/field_manager_r8.fh \
	include/fm_util.inc \
	include/fm_util_r4.fh \
	include/fm_util_r8.fh

field_manager_mod.$(FC_MODEXT): parse.inc fm_yaml_mod.$(FC_MODEXT) \
	include/field_manager.inc include/field_manager_r4.fh include/field_manager_r8.fh
fm_util_mod.$(FC_MODEXT): field_manager_mod.$(FC_MODEXT) \
	include/fm_util.inc include/fm_util_r4.fh include/fm_util_r8.fh

# Mod files are built and then installed as headers.
MODFILES = \
	field_manager_mod.$(FC_MODEXT) \
	fm_util_mod.$(FC_MODEXT) \
	fm_yaml_mod.$(FC_MODEXT)
BUILT_SOURCES = $(MODFILES)
nodist_include_HEADERS = $(MODFILES)

include $(top_srcdir)/mkmods.mk
