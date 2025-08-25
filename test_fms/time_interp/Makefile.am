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

# This is an automake file for the test_fms/time_interp directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_time_interp_r4 test_time_interp_r8 test_time_interp_external_r4 test_time_interp_external_r8

# These are the sources for the tests.
test_time_interp_r4_SOURCES = test_time_interp.F90
test_time_interp_r8_SOURCES = test_time_interp.F90
test_time_interp_external_r4_SOURCES = test_time_interp_external.F90
test_time_interp_external_r8_SOURCES = test_time_interp_external.F90

# filter out any added precision flags
# adds r8 flag, otherwise no-flag default is 4
test_time_interp_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -I$(MODDIR)
test_time_interp_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)
test_time_interp_external_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -I$(MODDIR)
test_time_interp_external_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -I$(MODDIR)

# Run the test programs.
TESTS = test_time_interp2.sh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh
# These files will also be included in the distribution.
EXTRA_DIST = test_time_interp2.sh

# Clean up
CLEANFILES = INPUT/* input.nml *.out* *.dpi *.spi *.dyn *.spl
