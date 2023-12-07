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
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = \
  test_data_override_init \
  test_get_grid_v1_r4 \
  test_get_grid_v1_r8 \
  test_data_override_r4 \
  test_data_override_r8 \
  test_data_override_ongrid_r4 \
  test_data_override_ongrid_r8

# This is the source code for the test.
test_data_override_init_SOURCES = test_data_override_init.F90
test_data_override_r4_SOURCES = test_data_override.F90
test_data_override_r8_SOURCES = test_data_override.F90

test_data_override_ongrid_r4_SOURCES = test_data_override_ongrid.F90
test_data_override_ongrid_r8_SOURCES = test_data_override_ongrid.F90

test_get_grid_v1_r4_SOURCES = test_get_grid_v1.F90
test_get_grid_v1_r8_SOURCES = test_get_grid_v1.F90

test_data_override_r4_CPPFLAGS = $(AM_CPPFLAGS) -DDO_TEST_KIND_=r4_kind
test_data_override_r8_CPPFLAGS = $(AM_CPPFLAGS) -DDO_TEST_KIND_=r8_kind

test_data_override_ongrid_r4_CPPFLAGS = $(AM_CPPFLAGS) -DDO_TEST_KIND_=r4_kind
test_data_override_ongrid_r8_CPPFLAGS = $(AM_CPPFLAGS) -DDO_TEST_KIND_=r8_kind

test_get_grid_v1_r4_CPPFLAGS = $(AM_CPPFLAGS) -DDO_TEST_KIND_=r4_kind
test_get_grid_v1_r8_CPPFLAGS = $(AM_CPPFLAGS) -DDO_TEST_KIND_=r8_kind

if SKIP_PARSER_TESTS
skipflag="skip"
else
skipflag=""
endif

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh
TESTS_ENVIRONMENT= test_input_path="@TEST_INPUT_PATH@" \
                   parser_skip=${skipflag}

# Run the test program.
TESTS = test_data_override2.sh test_data_override_init.sh

# Include these files with the distribution.
EXTRA_DIST = test_data_override2.sh test_data_override_init.sh

# Clean up
CLEANFILES = input.nml *.nc* *.out diag_table data_table data_table.yaml INPUT/* *.dpi *.spi *.dyn *.spl *-files/*
