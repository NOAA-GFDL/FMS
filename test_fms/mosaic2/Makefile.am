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
AM_CPPFLAGS = -I$(abs_top_srcdir)/test_fms/mosaic2 -I$(top_srcdir)/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_mosaic2_r4 test_mosaic2_r8 test_grid2_r4 test_grid2_r8

# This is the source code for the test
test_mosaic2_r4_SOURCES = write_files.F90 write_files.inc test_mosaic2.F90
test_grid2_r4_SOURCES   = write_files.F90 write_files.inc test_mosaic2.F90

test_mosaic2_r8_SOURCES = test_mosaic2.F90 write_files.inc write_files.F90
test_grid2_r8_SOURCES   = test_grid2.F90 write_files.inc write_files.F90

test_mosaic2_r4_CPPFLAGS=-DTEST_FMS_KIND_=4 -DWRITE_FILES_MOD_=write_files_r4 $(AM_CPPFLAGS)
test_grid2_r4_CPPFLAGS  =-DTEST_FMS_KIND_=4 -DWRITE_FILES_MOD_=write_files_r4 $(AM_CPPFLAGS)

test_mosaic2_r8_CPPFLAGS=-DTEST_FMS_KIND_=8 -DWRITE_FILES_MOD_=write_files_r8 $(AM_CPPFLAGS)
test_grid2_r8_CPPFLAGS  =-DTEST_FMS_KIND_=8 -DWRITE_FILES_MOD_=write_files_r8 $(AM_CPPFLAGS)

# compile helper module as well
%.o: %.mod
write_files_r8.mod: write_files.$(OBJEXT)
write_files_r4.mod: write_files.$(OBJEXT)
test_mosaic2_r4.$(OBJEXT): write_files_r4.mod
test_mosaic2_r8.$(OBJEXT): write_files_r8.mod
test_grid2_r4.$(OBJEXT): write_files_r4.mod
test_grid2_r8.$(OBJEXT): write_files_r8.mod

# These files are also included in the distribution.
EXTRA_DIST = test_mosaic2.sh

# Run the test program.
TESTS = test_mosaic2.sh

if SKIP_MOSAIC_TESTS
  TESTS_ENVIRONMENT = parser_skip=${skipflag}
endif

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
									$(abs_top_srcdir)/test_fms/tap-driver.sh
# Clean up
CLEANFILES = input.nml *.nc *.out *.dpi *.spi *.dyn *.spl *.mod
