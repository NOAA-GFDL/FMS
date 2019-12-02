# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez, Ed Hartnett

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Find the needed mod files.
AM_CPPFLAGS = -I${top_builddir}/horiz_interp -I${top_builddir}/fms \
-I${top_builddir}/constants -I${top_builddir}/mpp \
-I${top_builddir}/mosaic -I${top_builddir}/axis_utils

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_horiz_interp test2_horiz_interp

# These are the sources for the tests.
test_horiz_interp_SOURCES = test_horiz_interp.F90
test2_horiz_interp_SOURCES = test2_horiz_interp.F90

# Run the test programs.
TESTS = test_horiz_interp.sh test2_horiz_interp.sh

# These files will also be included in the distribution.
EXTRA_DIST = test_horiz_interp.sh test_horiz_interp.bats \
test2_horiz_interp.sh test2_horiz_interp.bats input_base.nml

# Clean up
CLEANFILES = input.nml *.out*
