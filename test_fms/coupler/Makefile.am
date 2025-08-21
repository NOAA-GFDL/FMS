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
