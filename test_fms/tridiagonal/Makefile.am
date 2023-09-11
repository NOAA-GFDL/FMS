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
check_PROGRAMS = test_tridiagonal4 test_tridiagonal8

# compiles test file with both kind sizes via macro
test_tridiagonal4.o:
	$(FC) $(AM_FCFLAGS) $(AM_CPPFLAGS) -DTEST_TRIDIAG_KIND=4 \
	    	-c $(top_srcdir)/test_fms/tridiagonal/test_tridiagonal.F90 \
			-o $@ $(AM_LDFLAGS) 

test_tridiagonal8.o:
	$(FC) $(AM_FCFLAGS) $(AM_CPPFLAGS) -DTEST_TRIDIAG_KIND=8 \
	    	-c $(top_srcdir)/test_fms/tridiagonal/test_tridiagonal.F90 \
			-o $@ $(AM_LDFLAGS) 

# Run the test program.
TESTS = test_tridiagonal.sh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# These files will be included in the distribution.
EXTRA_DIST = test_tridiagonal.sh

# Clean up
CLEANFILES = *.nml *.out* *.dpi *.spi *.dyn *.spl *.o test_tridiagonal4 test_tridiagonal8 test_tridiagonal
