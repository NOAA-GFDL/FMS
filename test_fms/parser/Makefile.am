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

# uramirez, Ed Hartnett

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I$(MODDIR)

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = parser_demo2 test_yaml_parser check_crashes parser_demo test_output_yaml

# This is the source code for the test.
test_yaml_parser_SOURCES = test_yaml_parser.F90
check_crashes_SOURCES = check_crashes.F90
parser_demo_SOURCES = parser_demo.F90
parser_demo2_SOURCES = parser_demo2.F90
test_output_yaml_SOURCES = test_output_yaml.F90

# Run the test program.
TESTS = test_yaml_parser.sh

# Include these files with the distribution.
EXTRA_DIST = test_yaml_parser.sh

# Define test file extensions and log driver
TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

if SKIP_PARSER_TESTS
skipflag="skip"
else
skipflag=""
endif

TESTS_ENVIRONMENT = parser_skip=${skipflag}
# Clean up
CLEANFILES = input.nml *.nc* *.out *.yaml *.dpi *.spi *.dyn *.spl
