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

# uramirez

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I$(MODDIR)

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_io_offloading test_mpp_redistribute test_metadata_transfer

# This is the source code for the test.
test_io_offloading_SOURCES = test_io_offloading.F90
test_mpp_redistribute_SOURCES = test_mpp_redistribute.F90
test_metadata_transfer_SOURCES= test_metadata_transfer.F90

# Run the test program.
TESTS = test_io_offloading.sh

# Define test file extensions and log driver
TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                   $(abs_top_srcdir)/test_fms/tap-driver.sh

# Include these files with the distribution.
EXTRA_DIST = test_io_offloading.sh

# Clean up
CLEANFILES = input.nml *.out
