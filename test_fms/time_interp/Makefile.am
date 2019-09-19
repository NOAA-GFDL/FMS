# This is an automake file for the test_fms/time_interp directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_time_interp_CPPFLAGS = -I${top_builddir}/time_manager
test_time_interp_CPPFLAGS += -I${top_builddir}/fms
test_time_interp_CPPFLAGS += -I${top_builddir}/time_interp

test_time_interp_external_CPPFLAGS = -I${top_builddir}/time_manager
test_time_interp_external_CPPFLAGS += -I${top_builddir}/fms
test_time_interp_external_CPPFLAGS += -I${top_builddir}/time_interp
test_time_interp_external_CPPFLAGS += -I${top_builddir}/mpp
test_time_interp_external_CPPFLAGS += -I${top_builddir}/axis_utils
test_time_interp_external_CPPFLAGS += -I${top_builddir}/horiz_interp
test_time_interp_external_CPPFLAGS += -I${top_builddir}/constants

# Copy over other needed files to the srcdir
EXTRA_DIST = test_time_interp.sh
EXTRA_DIST += test_time_interp.bats
EXTRA_DIST += test_time_interp_external.sh
EXTRA_DIST += test_time_interp_external.bats
EXTRA_DIST += input.nml

# Build this test program.
check_PROGRAMS = test_time_interp test_time_interp_external

# This is the source code for the test.
test_time_interp_SOURCES = test_time_interp.F90
test_time_interp_LDADD = ${top_builddir}/libFMS/libFMS.la

test_time_interp_external_SOURCES = test_time_interp_external.F90
test_time_interp_external_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_time_interp.sh test_time_interp_external.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml *.out*
