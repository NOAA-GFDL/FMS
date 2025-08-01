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

# This is an automake file for the test_fms/horiz_interp directory of
# the FMS package.

# uramirez, Ed Hartnett

# Find the needed mod files.
AM_CPPFLAGS = -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_horiz_interp_r4 test_horiz_interp_r8 test_create_xgrid_order2_r8

# These are the sources for the tests.
test_horiz_interp_r4_SOURCES = test_horiz_interp.F90
test_horiz_interp_r8_SOURCES = test_horiz_interp.F90
test_create_xgrid_order2_r8_SOURCES = test_create_xgrid_order2.F90

test_horiz_interp_r4_CPPFLAGS=-DHI_TEST_KIND_=4 -I$(MODDIR)
test_horiz_interp_r8_CPPFLAGS=-DHI_TEST_KIND_=8 -I$(MODDIR)
test_create_xgrid_order2_r8_CPPFLAGS=-DHI_TEST_KIND_=8 -I$(MODDIR)

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh
TESTS_ENVIRONMENT= test_input_path="@TEST_INPUT_PATH@"

# Run the test programs.
TESTS = test_horiz_interp2.sh test_create_xgrid_order2.sh

# These files will also be included in the distribution.
EXTRA_DIST = test_horiz_interp2.sh test_create_xgrid_order2.sh

# Clean up
CLEANFILES = input.nml *.out* *.dpi *.spi *.dyn *.spl
