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

# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I$(MODDIR)

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_string_utils

# This is the source code for the test.
test_string_utils_SOURCES = test_string_utils.F90

# Run the test program.
TESTS = test_string_utils.sh

# Define test file extensions and log driver
TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                   $(abs_top_srcdir)/test_fms/tap-driver.sh

# Include these files with the distribution.
EXTRA_DIST = test_string_utils.sh

# Clean up
CLEANFILES = input.nml *.out
