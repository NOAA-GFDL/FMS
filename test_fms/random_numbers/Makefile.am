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

# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# Find the fms and mpp mod files.
AM_CPPFLAGS = -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = \
    test_random_numbers_r4 \
    test_random_numbers_r8

# This is the source code for the test.
test_random_numbers_r4_SOURCES = test_random_numbers.F90
test_random_numbers_r8_SOURCES = test_random_numbers.F90

test_random_numbers_r4_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_RN_KIND_=r4_kind
test_random_numbers_r8_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_RN_KIND_=r8_kind

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# Run the test program.
TESTS = test_random_numbers.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = test_random_numbers.sh

# Clean up
CLEANFILES = input.nml *.out* *.nc* *.dpi *.spi *.dyn *.spl
