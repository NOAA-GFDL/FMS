# This is an automake file for the test_fms/time_interp directory of the FMS
# package.

# Find the fms_mod.mod file.
test_time_interp_CPPFLAGS = -I${top_srcdir}/time_manager
test_time_interp_CPPFLAGS += -I${top_srcdir}/fms
test_time_interp_CPPFLAGS += -I${top_srcdir}/time_interp

test_time_interp_external_CPPFLAGS = -I${top_srcdir}/time_manager
test_time_interp_external_CPPFLAGS += -I${top_srcdir}/fms
test_time_interp_external_CPPFLAGS += -I${top_srcdir}/time_interp
test_time_interp_external_CPPFLAGS += -I${top_srcdir}/mpp
test_time_interp_external_CPPFLAGS += -I${top_srcdir}/axis_utils
test_time_interp_external_CPPFLAGS += -I${top_srcdir}/horiz_interp
test_time_interp_external_CPPFLAGS += -I${top_srcdir}/constants

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
