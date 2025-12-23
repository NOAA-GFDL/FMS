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
