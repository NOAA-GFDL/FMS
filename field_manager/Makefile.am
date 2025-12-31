#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
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
	include/parse.inc \
	include/field_manager.inc \
	include/field_manager_r4.fh \
	include/field_manager_r8.fh \
	include/fm_util.inc \
	include/fm_util_r4.fh \
	include/fm_util_r8.fh

field_manager_mod.$(FC_MODEXT): include/parse.inc fm_yaml_mod.$(FC_MODEXT) \
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
