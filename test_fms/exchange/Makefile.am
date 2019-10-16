# This is an automake file for the test_fms/exchange directory of the FMS
# package.

# uramirez, Ed Hartnett

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)" builddir="$(top_builddir)"

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I${top_builddir}/mpp \
-I${top_builddir}/fms -I${top_builddir}/constants \
-I${top_builddir}/exchange -I${top_builddir}/mosaic \
-I${top_builddir}/coupler

# Link to FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_xgrid

# This is the source code for the test.
test_xgrid_SOURCES = test_xgrid.F90

# Run the test program.
TESTS = test_xgrid.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = test_xgrid.sh input.nml test_xgrid.bats

# Clean up
CLEANFILES = input.nml *.nc *.out
