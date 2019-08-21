# This is an automake file for the test_fms/interpolator directory of the FMS
# package.

# Find the fms_mod.mod file.
test_interpolator_CPPFLAGS = -I${top_srcdir}/mpp
test_interpolator_CPPFLAGS += -I${top_srcdir}/fms
test_interpolator_CPPFLAGS += -I${top_srcdir}/time_manager
test_interpolator_CPPFLAGS += -I${top_srcdir}/diag_manager
test_interpolator_CPPFLAGS += -I${top_srcdir}/interpolator
test_interpolator_CPPFLAGS += -I${top_srcdir}/constants
test_interpolator_CPPFLAGS += -I${top_srcdir}/time_interp

# Copy over other needed files to the srcdir
EXTRA_DIST = test_interpolator.sh
EXTRA_DIST +=  test_interpolator.bats
EXTRA_DIST += input.nml
EXTRA_DIST += diag_table

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
