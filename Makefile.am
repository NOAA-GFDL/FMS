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

# Directory to place all .mod files.  Convienient place to look for these
# files.
#include $(top_srcdir)/mkmods.mk

# Does the user want to build documentation?
if BUILD_DOCS
DOCS = docs
endif

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = \
  platform \
  mpp \
  tridiagonal \
  constants \
  constants4 \
  memutils \
  string_utils \
  fms2_io \
  mosaic2 \
  fms \
  parser \
  string_utils \
  affinity \
  grid_utils \
  time_manager \
  axis_utils \
  diag_manager \
  drifters \
  horiz_interp \
  time_interp \
  column_diagnostics \
  block_control \
  data_override \
  astronomy \
  field_manager \
  tracer_manager \
  coupler \
  diag_integral \
  monin_obukhov \
  interpolator \
  amip_interp \
  exchange \
  topography \
  sat_vapor_pres \
  random_numbers \
  . \
  libFMS \
  test_fms \
  ${DOCS}

include_HEADERS = include/file_version.h include/fms_platform.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = FMS.pc

## Build libFMS module
AM_CPPFLAGS = -I${top_srcdir}/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

noinst_LTLIBRARIES = libFMS_mod.la
libFMS_mod_la_SOURCES = libFMS.F90

fms.$(FC_MODEXT): .mods/*_mod.$(FC_MODEXT)

nodist_include_HEADERS = .mods/fms.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk

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


if COV
check-code-coverage: check
	-cd ${top_builddir}/test_fms && chmod +x intel_coverage.sh && ./intel_coverage.sh

.PHONY: check-code-coverage

clean-local:
	-rm -rf .mods coverage-data coverage-report test.nc
else
clean-local:
	-rm -rf .mods test.nc
endif

install-data-hook:
	@echo ''
	@echo '+-------------------------------------------------------------+'
	@echo '| Congratulations! You have successfully installed the FMS    |'
	@echo '| Fortran library.                                            |'
	@echo '|                                                             |'
	@echo '| After the installed include and link paths have been        |'
	@echo '| specified, code using FMS should be compiled using the      |'
	@echo '| "-lFMS" flag.                                               |'
	@echo '|                                                             |'
	@echo '| FMS is developed and maintained at the GFDL publicly on     |'
	@echo '| Github. To report an issue or view available documentation, |'
	@echo '| please see our page: https://www.github.com/NOAA-GFDL/FMS   |'
	@echo '+-------------------------------------------------------------+'
	@echo ''
