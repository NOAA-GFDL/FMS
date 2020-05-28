# This is the main automake file for FMS.
# Ed Hartnett 2/21/2019

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = include platform constants tridiagonal mpp memutils fms \
affinity mosaic time_manager diag_manager drifters axis_utils \
horiz_interp time_interp column_diagnostics block_control \
data_override astronomy field_manager coupler diag_integral \
monin_obukhov interpolator fft amip_interp oda_tools exchange \
topography tracer_manager station_data sat_vapor_pres random_numbers \
libFMS test_fms

# Prepare CMake files for installation.  This is to help
# packages build using CMake to more easily use the libFMS
# library.
#
# A simple edit command
edit = sed \
         -e 's|@PACKAGE_VERSION[@]|$(pkgdatadir)|g' \
         -e 's|@PACKAGE_CMAKE_INSTALL_PREFIX[@]|$(prefix)|g' \
         -e 's|@PACKAGE_CMAKE_INSTALL_INCLUDEDIR[@]|$(includedir)|g' \
         -e 's|@PACKAGE_CMAKE_INSTALL_LIBDIR[@]|$(libdir)|g' \
         -e '/@PACKAGE_INIT[@]/ {' -e 'r $(top_srcdir)/_package_init_cmake' -e  'd' -e '}'
FMSConfig.cmake FMSConfigVersion.cmake: Makefile
	rm -f $@ $@.tmp
	$(edit) $@.in >$@.tmp
	mv $@.tmp $@
FMSConfig.cmake: FMSConfig.cmake.in
FMSConfigVersion.cmake: FMSConfigVersion.cmake.in

cmakedir = $(libdir)
cmake_DATA = FMSConfig.cmake \
  FMSConfigVersion.cmake
 
EXTRA_DIST = README.md _package_init_cmake

CLEANFILES = *__genmod.mod *__genmod.f90
