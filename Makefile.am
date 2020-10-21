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

# This is the main automake file for FMS.
# Ed Hartnett 2/21/2019

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# Does the user want to build documentation?
if BUILD_DOCS
DOCS = docs
endif

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = include platform constants tridiagonal mpp memutils fms2_io \
fms affinity mosaic time_manager axis_utils diag_manager drifters \
horiz_interp time_interp column_diagnostics block_control \
data_override astronomy field_manager coupler diag_integral \
monin_obukhov interpolator amip_interp exchange topography \
tracer_manager sat_vapor_pres random_numbers libFMS test_fms ${DOCS}

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = FMS.pc

# Prepare CMake files for installation.  This is to help
# packages build using CMake to more easily use the libFMS
# library.
#
# A simple edit command
#edit = sed \
#         -e 's|@PACKAGE_VERSION[@]|$(pkgdatadir)|g' \
#         -e 's|@PACKAGE_CMAKE_INSTALL_PREFIX[@]|$(prefix)|g' \
#         -e 's|@PACKAGE_CMAKE_INSTALL_INCLUDEDIR[@]|$(includedir)|g' \
#         -e 's|@PACKAGE_CMAKE_INSTALL_LIBDIR[@]|$(libdir)|g' \
#         -e '/@PACKAGE_INIT[@]/ {' -e 'r $(top_srcdir)/_package_init_cmake' -e  'd' -e '}'
#FMSConfig.cmake FMSConfigVersion.cmake: Makefile
#	rm -f $@ $@.tmp
#	$(edit) $@.in >$@.tmp
#	mv $@.tmp $@
#FMSConfig.cmake: FMSConfig.cmake.in
#FMSConfigVersion.cmake: FMSConfigVersion.cmake.in
#
#cmakedir = $(libdir)
#cmake_DATA = FMSConfig.cmake \
#  FMSConfigVersion.cmake

#clean-local:
#	-rm -f FMSConfig.cmake FMSConfigVersion.cmake
EXTRA_DIST = README.md #_package_init_cmake

CLEANFILES = *__genmod.mod *__genmod.f90
