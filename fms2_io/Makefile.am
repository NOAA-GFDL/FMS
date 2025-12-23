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

# This is an automake file for the fms2_io directory of the FMS package.

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/fms2_io/include -I$(top_srcdir)/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfms2_io.la

# Each convenience library depends on its source.
libfms2_io_la_SOURCES = \
  fms2_io.F90 \
  fms_io_utils.F90 \
  netcdf_io.F90 \
  fms_netcdf_domain_io.F90 \
  fms_netcdf_unstructured_domain_io.F90 \
  blackboxio.F90 \
  include/array_utils_char.inc \
  include/compressed_read.inc \
  include/compute_global_checksum.inc \
  include/domain_write.inc \
  include/get_data_type_string.inc \
  include/get_variable_attribute.inc \
  include/netcdf_read_data.inc \
  include/register_domain_restart_variable.inc \
  include/register_unstructured_domain_restart_variable.inc \
  include/unstructured_domain_read.inc \
  include/array_utils.inc \
  include/compressed_write.inc \
  include/domain_read.inc \
  include/get_global_attribute.inc \
  include/netcdf_add_restart_variable.inc \
  include/netcdf_write_data.inc \
  include/register_global_attribute.inc \
  include/register_variable_attribute.inc \
  include/unstructured_domain_write.inc \
  include/gather_data_bc.inc \
  include/scatter_data_bc.inc \
  include/unpack_data.inc

# Some mods are dependant on other mods in this dir.
fms2_io_mod.$(FC_MODEXT): fms_io_utils_mod.$(FC_MODEXT) netcdf_io_mod.$(FC_MODEXT) fms_netcdf_domain_io_mod.$(FC_MODEXT) \
                 fms_netcdf_unstructured_domain_io_mod.$(FC_MODEXT) blackboxio.$(FC_MODEXT)
fms_io_utils_mod.$(FC_MODEXT): include/array_utils.inc include/array_utils_char.inc \
                      include/get_data_type_string.inc
netcdf_io_mod.$(FC_MODEXT): fms_io_utils_mod.$(FC_MODEXT) include/netcdf_add_restart_variable.inc include/netcdf_read_data.inc \
                   include/netcdf_write_data.inc include/register_global_attribute.inc \
                   include/register_variable_attribute.inc include/get_global_attribute.inc \
                   include/get_variable_attribute.inc include/compressed_write.inc include/compressed_read.inc \
                   include/gather_data_bc.inc include/scatter_data_bc.inc include/unpack_data.inc
fms_netcdf_domain_io_mod.$(FC_MODEXT): fms_io_utils_mod.$(FC_MODEXT) netcdf_io_mod.$(FC_MODEXT) include/register_domain_restart_variable.inc \
                              include/domain_read.inc include/domain_write.inc include/compute_global_checksum.inc
fms_netcdf_unstructured_domain_io_mod.$(FC_MODEXT): fms_io_utils_mod.$(FC_MODEXT) netcdf_io_mod.$(FC_MODEXT) \
                                           include/register_unstructured_domain_restart_variable.inc \
                                           include/unstructured_domain_read.inc include/unstructured_domain_write.inc
blackboxio.$(FC_MODEXT): fms_io_utils_mod.$(FC_MODEXT) netcdf_io_mod.$(FC_MODEXT) fms_netcdf_domain_io_mod.$(FC_MODEXT) \
                fms_netcdf_unstructured_domain_io_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
MODFILES = \
  fms2_io_mod.$(FC_MODEXT) \
  fms_io_utils_mod.$(FC_MODEXT) \
  netcdf_io_mod.$(FC_MODEXT) \
  fms_netcdf_domain_io_mod.$(FC_MODEXT) \
  fms_netcdf_unstructured_domain_io_mod.$(FC_MODEXT) \
  blackboxio.$(FC_MODEXT)
BUILT_SOURCES = $(MODFILES)
nodist_include_HEADERS = $(MODFILES)

include $(top_srcdir)/mkmods.mk
