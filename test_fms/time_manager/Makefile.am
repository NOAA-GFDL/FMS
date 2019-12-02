# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez, Ed Hartnett

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/mpp -I${top_builddir}/fms \
-I${top_builddir}/constants -I${top_builddir}/time_manager

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_time_manager

# This is the source code for the test.
test_time_manager_SOURCES = test_time_manager.F90

# Run the test program.
TESTS = test_time_manager.sh

# These files will be included in the distribution.
EXTRA_DIST = test_time_manager.sh test_time_manager.bats \
input_base.nml

# Clean up
CLEANFILES = input.nml *.out*
