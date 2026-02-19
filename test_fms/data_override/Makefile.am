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

# uramirez, Ed Hartnett

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I$(MODDIR) -I$(top_srcdir)/include -I$(top_srcdir)/test_fms/data_override/include -I$(top_builddir)/test_fms/common

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la $(top_builddir)/test_fms/common/libtest_fms.a

# Build this test program.
check_PROGRAMS = \
  test_data_override_init \
  test_get_grid_v1_r4 \
  test_get_grid_v1_r8 \
  test_data_override_r4 \
  test_data_override_r8 \
  test_data_override_ongrid

# This is the source code for the test.
test_data_override_init_SOURCES = test_data_override_init.F90
test_data_override_r4_SOURCES = test_data_override.F90
test_data_override_r8_SOURCES = test_data_override.F90

test_data_override_ongrid_SOURCES = test_data_override_ongrid.F90 \
				    include/test_data_override_ongrid.inc \
				    include/test_data_override_ongrid_r4.fh \
				    include/test_data_override_ongrid_r8.fh

test_get_grid_v1_r4_SOURCES = test_get_grid_v1.F90
test_get_grid_v1_r8_SOURCES = test_get_grid_v1.F90

test_data_override_r4_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_FMS_KIND_=r4_kind
test_data_override_r8_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_FMS_KIND_=r8_kind

test_get_grid_v1_r4_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_FMS_KIND_=r4_kind
test_get_grid_v1_r8_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_FMS_KIND_=r8_kind


if USING_YAML
  TESTS_ENVIRONMENT = parser_skip=""
else
  TESTS_ENVIRONMENT = parser_skip=skip
endif

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh
TESTS_ENVIRONMENT+= test_input_path="@TEST_INPUT_PATH@"

# Run the test program.

TESTS = test_data_override2.sh test_data_override_init.sh test_data_override2_mono.sh test_data_override2_ongrid.sh \
        test_data_override2_scalar.sh test_data_override_weights.sh test_data_override_ensembles.sh \
        test_data_override_multi_file.sh

# Include these files with the distribution.
EXTRA_DIST = test_data_override2.sh test_data_override_init.sh test_data_override2_mono.sh test_data_override2_ongrid.sh \
             test_data_override2_scalar.sh test_data_override_weights.sh test_data_override_ensembles.sh \
             test_data_override_multi_file.sh

# Clean up
CLEANFILES = input.nml *.nc* *.out diag_table data_table data_table.yaml INPUT/* *.dpi *.spi *.dyn *.spl *-files/*
