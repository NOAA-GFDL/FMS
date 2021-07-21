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

# This is an automake file for the test_fms/mpp_io directory of the
# FMS package.

# uramirez, Ed Hartnett, Ryan Mulhall

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_mpp_io \
	test_io_R4_R8              \
	test_io_mosaic_R4_R8

# This is the source code for the test.
test_mpp_io_SOURCES = test_mpp_io.F90
test_io_R4_R8_SOURCES = test_io_R4_R8.F90
test_io_mosaic_R4_R8_SOURCES = test_io_mosaic_R4_R8.F90

# Run the test program.
TESTS = test_mpp_io2.sh \
	test_io_R4_R8.sh      \
	test_io_mosaic_R4_R8.sh

# These files will also be distributed.
EXTRA_DIST = test_mpp_io2.sh test_io_R4_R8.sh test_io_mosaic_R4_R8.sh input_base.nml

# Clean up
CLEANFILES = input.nml *.nc* *.out
