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
AM_CPPFLAGS = -I${top_srcdir}/include -I$(MODDIR)

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_coupler_2d_r4 test_coupler_2d_r8 \
                 test_coupler_3d_r4 test_coupler_3d_r8 \
                 test_coupler_types_r4 test_coupler_types_r8 \
                 test_atmos_ocean_fluxes_r4 test_atmos_ocean_fluxes_r8 \
                 test_gex

# This is the source code for the test.
test_coupler_2d_r4_SOURCES = test_coupler_2d.F90 test_coupler_utils.inc
test_coupler_2d_r8_SOURCES = test_coupler_2d.F90 test_coupler_utils.inc
test_coupler_3d_r4_SOURCES = test_coupler_3d.F90 test_coupler_utils.inc
test_coupler_3d_r8_SOURCES = test_coupler_3d.F90 test_coupler_utils.inc
test_coupler_types_r4_SOURCES = test_coupler_types.F90 test_coupler_utils.inc
test_coupler_types_r8_SOURCES = test_coupler_types.F90 test_coupler_utils.inc
test_atmos_ocean_fluxes_r4_SOURCES = test_atmos_ocean_fluxes.F90
test_atmos_ocean_fluxes_r8_SOURCES = test_atmos_ocean_fluxes.F90
test_gex_SOURCES = test_gex.F90

# add macro flags for kind and r4/r8 child type
test_coupler_2d_r4_CPPFLAGS=-DTEST_FMS_KIND_=r4_kind -DFMS_TEST_BC_TYPE_=bc_r4 -I$(MODDIR)
test_coupler_2d_r8_CPPFLAGS=-DTEST_FMS_KIND_=r8_kind -DFMS_TEST_BC_TYPE_=bc -I$(MODDIR)
test_coupler_3d_r4_CPPFLAGS=-DTEST_FMS_KIND_=r4_kind -DFMS_TEST_BC_TYPE_=bc_r4 -I$(MODDIR)
test_coupler_3d_r8_CPPFLAGS=-DTEST_FMS_KIND_=r8_kind -DFMS_TEST_BC_TYPE_=bc -I$(MODDIR)
test_coupler_types_r4_CPPFLAGS=-DTEST_FMS_KIND_=r4_kind -DFMS_TEST_BC_TYPE_=bc_r4 -I$(MODDIR)
test_coupler_types_r8_CPPFLAGS=-DTEST_FMS_KIND_=r8_kind -DFMS_TEST_BC_TYPE_=bc -I$(MODDIR)
test_atmos_ocean_fluxes_r4_CPPFLAGS=-DTEST_FMS_KIND_=r4_kind -DFMS_TEST_BC_TYPE_=bc_r4 -I$(MODDIR)
test_atmos_ocean_fluxes_r8_CPPFLAGS=-DTEST_FMS_KIND_=r8_kind -DFMS_TEST_BC_TYPE_=bc -I$(MODDIR)
test_gex_CPPFLAGS = -I$(MODDIR)

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# Run the test program.
TESTS = test_coupler.sh test_gex.sh

if USING_YAML
  TESTS_ENVIRONMENT = parser_skip=""
else
  TESTS_ENVIRONMENT = parser_skip=skip
endif

# Include these files with the distribution.
EXTRA_DIST = test_coupler.sh test_gex.sh

# Clean up
CLEANFILES = input.nml *.nc* *.out *.dpi *.spi *.dyn *.spl *_table* INPUT/*.nc test_coupler-files/* test_gex-files/* test_coupler-files/*/* test_gex-files/*/*
