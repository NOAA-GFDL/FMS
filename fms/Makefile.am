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

# This is an automake file for the fms directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfms.la

# Each convenience library depends on its source.
libfms_la_SOURCES = \
  fms.F90 \
  fms_io.F90 \
  fms_io_unstructured_field_exist.inc \
  fms_io_unstructured_get_file_name.inc \
  fms_io_unstructured_register_restart_axis.inc \
  fms_io_unstructured_setup_one_field.inc read_data_4d.inc \
  fms_io_unstructured_file_unit.inc \
  fms_io_unstructured_get_file_unit.inc \
  fms_io_unstructured_register_restart_field.inc \
  read_data_2d.inc \
  write_data.inc \
  fms_io_unstructured_get_field_size.inc \
  fms_io_unstructured_read.inc \
  fms_io_unstructured_save_restart.inc \
  read_data_3d.inc

fms_mod.$(FC_MODEXT): fms_io_mod.$(FC_MODEXT)
fms_io_mod.$(FC_MODEXT): fms_io_unstructured_field_exist.inc \
  fms_io_unstructured_get_file_name.inc \
  fms_io_unstructured_register_restart_axis.inc \
  fms_io_unstructured_setup_one_field.inc read_data_4d.inc \
  fms_io_unstructured_file_unit.inc \
  fms_io_unstructured_get_file_unit.inc \
  fms_io_unstructured_register_restart_field.inc \
  read_data_2d.inc \
  write_data.inc \
  fms_io_unstructured_get_field_size.inc \
  fms_io_unstructured_read.inc \
  fms_io_unstructured_save_restart.inc \
  read_data_3d.inc

# Mod files are built and then installed as headers.
MODFILES = \
  fms_io_mod.$(FC_MODEXT) \
  fms_mod.$(FC_MODEXT)
BUILT_SOURCES = $(MODFILES)
nodist_include_HEADERS = $(MODFILES) $(FMS_INC_FILES)

include $(top_srcdir)/mkmods.mk
