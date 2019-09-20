# This is an automake file for the test_fms/data_override directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/time_manager

# Copy over other needed files to the srcdir
EXTRA_DIST =test_time_manager.sh
EXTRA_DIST += test_time_manager.bats
EXTRA_DIST += input_base.nml

# Build this test program.
check_PROGRAMS = test_time_manager

# This is the source code for the test.
test_time_manager_SOURCES = test_time_manager.F90
test_time_manager_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_time_manager.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml *.out*
