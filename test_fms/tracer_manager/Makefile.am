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

# This is an automake file for the test_fms/tracer_manager directory of
# the FMS package.

# uramirez, Ed Hartnett

# Find the needed mod files.
AM_CPPFLAGS = -I$(MODDIR) -I$(top_srcdir)/include

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_tracer_manager_r4 test_tracer_manager_r8

# This is the source code for the test.
test_tracer_manager_r4_SOURCES = test_tracer_manager.F90
test_tracer_manager_r8_SOURCES = test_tracer_manager.F90

test_tracer_manager_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -I$(MODDIR)
test_tracer_manager_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)

if USING_YAML
  TESTS_ENVIRONMENT = parser_skip=""
else
  TESTS_ENVIRONMENT = parser_skip=skip
endif

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) $(abs_top_srcdir)/test_fms/tap-driver.sh

# Run the test program.
TESTS = test_tracer_manager2.sh

# These files will also be included in the distribution.
EXTRA_DIST = test_tracer_manager2.sh

# Clean up
CLEANFILES = input.nml *.out* field_table *.dpi *.spi *.dyn *.spl *.yaml
