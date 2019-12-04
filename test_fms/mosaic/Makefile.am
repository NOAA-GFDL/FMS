# This is an automake file for the test_fms/mosaic directory of the FMS
# package.

# uramirez, Ed Hartnett

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Find the needed mod and include files.
AM_CPPFLAGS = -I${top_builddir}/mosaic -I${top_builddir}/fms \
-I${top_srcdir}/include -I${top_builddir}/mpp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_mosaic

# This is the source code for the test.
test_mosaic_SOURCES = test_mosaic.F90

# Run the test program.
TESTS = test_mosaic.sh

# These files are also included in the distribution.
EXTRA_DIST = test_mosaic.sh test_mosaic.bats input.nml
