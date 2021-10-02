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

# This is an automake file for the test_fms/fms2_io directory of the FMS
# package.

# J.Liptak, R.Menzel, U.Ramirez

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_get_is_valid test_file_appendix test_fms2_io test_atmosphere_io test_io_simple test_io_with_mask test_global_att \
                 test_bc_restart test_get_mosaic_tile_grid test_read_ascii_file

# This is the source code for the test.
test_get_is_valid_SOURCES = test_get_is_valid.F90
test_fms2_io_SOURCES = argparse.F90 test_fms2_io.F90 create_atmosphere_domain.inc \
                       create_land_domain.inc create_ocean_domain.inc \
                       ocean_restart_file_test.inc land_compressed_restart_file_test.inc \
                       land_unstructured_restart_file_test.inc
test_atmosphere_io_SOURCES = argparse.F90 setup.F90 test_atmosphere_io.F90 \
                             atmosphere_restart_file_test.inc
test_io_simple_SOURCES = test_io_simple.F90 argparse.F90 setup.F90
test_global_att_SOURCES = test_global_att.F90
test_io_with_mask_SOURCES=test_io_with_mask.F90
test_bc_restart_SOURCES=test_bc_restart.F90
test_get_mosaic_tile_grid_SOURCES=test_get_mosaic_tile_grid.F90
test_read_ascii_file_SOURCES=test_read_ascii_file.F90
test_file_appendix_SOURCES=test_file_appendix.F90

EXTRA_DIST = test_bc_restart.sh test_fms2_io.sh test_atmosphere_io.sh test_io_simple.sh test_global_att.sh test_io_with_mask.sh test_read_ascii_file.sh

argparse.mod: argparse.$(OBJEXT)
setup.mod: setup.$(OBJEXT)
setup.$(OBJEXT): argparse.mod
test_atmosphere_io.$(OBJEXT): setup.mod
test_io_simple.$(OBJEXT): setup.mod
test_fms2_io.$(OBJEXT): argparse.mod

# Run the test program.
TESTS = test_bc_restart.sh test_fms2_io.sh test_atmosphere_io.sh test_io_simple.sh test_io_with_mask.sh test_global_att.sh test_read_ascii_file.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)";

CLEANFILES = *.mod *.nc *.nc.* input.nml logfile.000000.out the_mask ascii_test1
