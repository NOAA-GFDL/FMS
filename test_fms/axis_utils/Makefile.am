# This is an automake file for the test_fms/data_override directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/fms -I${top_builddir}/mpp
LDADD = ${top_builddir}/libFMS/libFMS.la

# Copy over other needed files to the srcdir
EXTRA_DIST = test_axis_utils.sh test_axis_utils.bats input.nml

# Build this test program.
check_PROGRAMS = test_axis_utils

# This is the source code for the test.
test_axis_utils_SOURCES = test_axis_utils.F90

# Run the test program.
TESTS = test_axis_utils.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml *.out*
