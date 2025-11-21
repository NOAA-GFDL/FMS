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

# This is an automake file for the test_fms/topography directory of the
# FMS package.

# Caitlyn McAllister

# Find the fms and mpp mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = \
  test_topography_r4 \
  test_topography_r8

# This is the source code for the test.
test_topography_r4_SOURCES = test_topography.F90
test_topography_r8_SOURCES = test_topography.F90

# Set r4_kind and r8_kind
test_topography_r4_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_FMS_KIND_=r4_kind
test_topography_r8_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_FMS_KIND_=r8_kind

# Run the test program.
TESTS = test_topography.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = test_topography.sh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# Clean up
CLEANFILES = input.nml *.out* *.dpi *.spi *.dyn *.spl *.nc
