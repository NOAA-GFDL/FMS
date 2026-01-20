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

# This is an automake file for the test_fms/time_interp directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_time_interp_r4 test_time_interp_r8 \
                 test_time_interp_external_r4 test_time_interp_external_r8 \
                 test_time_interp_conservative_hi_r4 test_time_interp_conservative_hi_r8

# These are the sources for the tests.
test_time_interp_r4_SOURCES = test_time_interp.F90
test_time_interp_r8_SOURCES = test_time_interp.F90
test_time_interp_external_r4_SOURCES = test_time_interp_external.F90
test_time_interp_external_r8_SOURCES = test_time_interp_external.F90
test_time_interp_conservative_hi_r4_SOURCES = test_time_interp_conservative_hi.F90
test_time_interp_conservative_hi_r8_SOURCES = test_time_interp_conservative_hi.F90

# filter out any added precision flags
# adds r8 flag, otherwise no-flag default is 4
test_time_interp_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -I$(MODDIR)
test_time_interp_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)
test_time_interp_external_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -I$(MODDIR)
test_time_interp_external_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)
test_time_interp_conservative_hi_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -I$(MODDIR)
test_time_interp_conservative_hi_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)

# Run the test programs.
TESTS = test_time_interp2.sh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh
# These files will also be included in the distribution.
EXTRA_DIST = test_time_interp2.sh

# Clean up
CLEANFILES = INPUT/* input.nml *.out* *.dpi *.spi *.dyn *.spl
