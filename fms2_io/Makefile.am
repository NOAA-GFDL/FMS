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

# This is an automake file for the fms2_io directory of the FMS package.

# Descend into include directory.
SUBDIRS = include

# These are the include files that contain fortran code.
INCFILES = include/array_utils_char.inc include/compressed_read.inc	\
include/compute_global_checksum.inc include/domain_write.inc		\
include/get_data_type_string.inc include/get_variable_attribute.inc	\
include/netcdf_read_data.inc						\
include/register_domain_restart_variable.inc				\
include/register_unstructured_domain_restart_variable.inc		\
include/unstructured_domain_read.inc include/array_utils.inc		\
include/compressed_write.inc include/domain_read.inc			\
include/get_checksum.inc include/get_global_attribute.inc		\
include/netcdf_add_restart_variable.inc include/netcdf_write_data.inc	\
include/register_global_attribute.inc					\
include/register_variable_attribute.inc					\
include/unstructured_domain_write.inc

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_srcdir}/fms2_io/include

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfms2_io.la libfms_io_utils.la libnetcdf_io.la	\
libfms_netcdf_domain_io.la libfms_netcdf_unstructured_domain_io.la	\
libblackboxio.la

# Each convenience library depends on its source.
libfms2_io_la_SOURCES = fms2_io.F90
libfms_io_utils_la_SOURCES = fms_io_utils.F90
libnetcdf_io_la_SOURCES = netcdf_io.F90
libfms_netcdf_domain_io_la_SOURCES = fms_netcdf_domain_io.F90
libfms_netcdf_unstructured_domain_io_la_SOURCES = fms_netcdf_unstructured_domain_io.F90
libblackboxio_la_SOURCES = blackboxio.F90

# Each mod file depends on the .o file.
fms2_io_mod.mod: fms2_io.lo
fms_io_utils_mod.mod: fms_io_utils.lo
netcdf_io_mod.mod: netcdf_io.lo
fms_netcdf_domain_io_mod.mod: fms_netcdf_domain_io.lo
fms_netcdf_unstructured_domain_io_mod.mod: fms_netcdf_unstructured_domain_io.lo
blackboxio.mod: blackboxio.lo

# Some mods are dependant on other mods in this dir.
fms2_io.lo: fms_io_utils_mod.mod netcdf_io_mod.mod fms_netcdf_domain_io_mod.mod \
fms_netcdf_unstructured_domain_io_mod.mod blackboxio.mod $(INCFILES)
netcdf_io.lo: fms_io_utils_mod.mod $(INCFILES)
fms_netcdf_domain_io.lo: fms_io_utils_mod.mod netcdf_io_mod.mod $(INCFILES)
fms_netcdf_unstructured_domain_io.lo: fms_io_utils_mod.mod netcdf_io_mod.mod $(INCFILES)
blackboxio.lo: fms_io_utils_mod.mod netcdf_io_mod.mod fms_netcdf_domain_io_mod.mod \
fms_netcdf_unstructured_domain_io_mod.mod $(INCFILES)

# Mod files are built and then installed as headers.
MODFILES = fms2_io_mod.mod fms_io_utils_mod.mod netcdf_io_mod.mod	\
fms_netcdf_domain_io_mod.mod						\
fms_netcdf_unstructured_domain_io_mod.mod blackboxio.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

CLEANFILES = *.mod
