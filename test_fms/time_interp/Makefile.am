# This is an automake file for the test_fms/time_interp directory of the FMS
# package.

# uramirez, Ed Hartnett

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/time_manager -I${top_builddir}/fms \
-I${top_builddir}/time_interp -I${top_builddir}/mpp \
-I${top_builddir}/axis_utils -I${top_builddir}/horiz_interp \
-I${top_builddir}/constants

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_time_interp test_time_interp_external

# These are the sources for the tests.
test_time_interp_SOURCES = test_time_interp.F90
test_time_interp_external_SOURCES = test_time_interp_external.F90

# Run the test programs.
TESTS = test_time_interp.sh test_time_interp_external.sh

# These files will also be included in the distribution.
EXTRA_DIST = test_time_interp.sh test_time_interp.bats \
test_time_interp_external.sh test_time_interp_external.bats \
input_base.nml

# Clean up
CLEANFILES = input.nml *.out*
