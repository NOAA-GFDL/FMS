# This is an automake file for the test_fms/interpolator directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_interpolator_CPPFLAGS = -I${top_builddir}/mpp
test_interpolator_CPPFLAGS += -I${top_builddir}/fms
test_interpolator_CPPFLAGS += -I${top_builddir}/time_manager
test_interpolator_CPPFLAGS += -I${top_builddir}/diag_manager
test_interpolator_CPPFLAGS += -I${top_builddir}/interpolator
test_interpolator_CPPFLAGS += -I${top_builddir}/constants
test_interpolator_CPPFLAGS += -I${top_builddir}/time_interp

# Copy over other needed files to the srcdir
EXTRA_DIST = test_interpolator.sh
EXTRA_DIST +=  test_interpolator.bats
EXTRA_DIST += input_base.nml
EXTRA_DIST += diag_table_base

# Build this test program.
check_PROGRAMS = test_interpolator

# This is the source code for the test.
test_interpolator_SOURCES = test_interpolator.F90
test_interpolator_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_interpolator.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml diag_table
