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
check_PROGRAMS = test_icebergs

# This is the source code for the test.
test_icebergs_SOURCES=test_icebergs.F90

EXTRA_DIST = test_bc_restart.sh test_fms2_io.sh test_atmosphere_io.sh test_io_simple.sh test_global_att.sh test_io_with_mask.sh test_read_ascii_file.sh

# fixes issue with make implicitly trying to use m2c for argparse
%.o: %.mod

# Run the test program.
#TESTS = test_bc_restart.sh test_fms2_io.sh test_atmosphere_io.sh test_io_simple.sh test_io_with_mask.sh test_global_att.sh test_read_ascii_file.sh

# Set srcdir as evironment variable to be reference in the job script
TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh
TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)"

CLEANFILES = *.mod *.nc *.nc.* input.nml logfile.000000.out the_mask ascii_test1 *.dpi *.spi *.dyn *.spl *-files/*
