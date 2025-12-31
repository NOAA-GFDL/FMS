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

# This is an automake file for the test_fms/horiz_interp directory of
# the FMS package.

# uramirez, Ed Hartnett

# Find the needed mod files.
AM_CPPFLAGS = -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_horiz_interp_r4 \
                 test_horiz_interp_r8 \
                 test_create_xgrid_order2_r8 \
                 test_horiz_interp_read \
                 test_horiz_interp_fregrid_compatibility

# These are the sources for the tests.
test_horiz_interp_r4_SOURCES = test_horiz_interp.F90
test_horiz_interp_r8_SOURCES = test_horiz_interp.F90
test_create_xgrid_order2_r8_SOURCES = test_create_xgrid_order2.F90
test_horiz_interp_read_SOURCES = test_horiz_interp_read.F90
test_horiz_interp_fregrid_compatibility_SOURCES = test_horiz_interp_fregrid_compatibility.F90

test_horiz_interp_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -I$(MODDIR)
test_horiz_interp_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)
test_create_xgrid_order2_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh
TESTS_ENVIRONMENT= test_input_path="@TEST_INPUT_PATH@"

# Run the test programs.
TESTS = test_horiz_interp2.sh test_create_xgrid_order2.sh

# These files will also be included in the distribution.
EXTRA_DIST = test_horiz_interp2.sh test_create_xgrid_order2.sh

# Clean up
CLEANFILES = input.nml *.out* *.dpi *.spi *.dyn *.spl *nc
