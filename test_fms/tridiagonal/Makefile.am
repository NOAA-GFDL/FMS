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

# This is an automake file for the test_fms/tridiagonal directory of the FMS
# package.

# Ryan Mulhall

# Find the .mod directory
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_tridiagonal_r4 test_tridiagonal_r8

# compiles test file with both kind sizes via macro
test_tridiagonal_r4_SOURCES=test_tridiagonal.F90
test_tridiagonal_r8_SOURCES=test_tridiagonal.F90

test_tridiagonal_r4_CPPFLAGS=-DTEST_FMS_KIND_=r4_kind -I$(MODDIR)
test_tridiagonal_r8_CPPFLAGS=-DTEST_FMS_KIND_=r8_kind -I$(MODDIR)

# Run the test program.
TESTS = test_tridiagonal.sh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# These files will be included in the distribution.
EXTRA_DIST = test_tridiagonal.sh

# Clean up
CLEANFILES = *.nml *.out* *.dpi *.spi *.dyn *.spl *.o test_tridiagonal4 test_tridiagonal8 test_tridiagonal
