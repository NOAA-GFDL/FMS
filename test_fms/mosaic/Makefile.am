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

# This is an automake file for the test_fms/mosaic directory of the
# FMS package.

# uramirez, Ed Hartnett

# Find the needed mod and include files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(MODDIR) -I./INPUT -I$(top_srcdir)/mosaic

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_mosaic test_grid2

# This is the source code for the test
test_mosaic_SOURCES = test_mosaic.F90 write_files.F90
test_grid2_SOURCES = test_grid2.F90 write_files.F90

# These files are also included in the distribution.
EXTRA_DIST = test_mosaic2.sh

%.o: %.mod
write_files.mod: write_files.$(OBJEXT)
test_mosaic.$(OBJEXT): write_files.mod
test_grid2.$(OBJEXT): write_files.mod

# Run the test program.
TESTS = test_mosaic2.sh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
									$(abs_top_srcdir)/test_fms/tap-driver.sh

# Clean up
CLEANFILES = input.nml *.nc *.out *.dpi *.spi *.dyn *.spl
