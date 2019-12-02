# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez, Ed Hartnett

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script.
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)" builddir="$(top_builddir)"

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I${top_builddir}/mpp \
-I${top_builddir}/fms -I${top_builddir}/constants \
-I${top_builddir}/time_manager -I${top_builddir}/diag_manager \
-I${top_builddir}/data_override

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_data_override

# This is the source code for the test.
test_data_override_SOURCES = test_data_override.F90

# Run the test program.
TESTS = test_data_override.sh

# Include these files with the distribution.
EXTRA_DIST = test_data_override.sh test_data_override.bats \
input.nml_base diag_table data_table

# Clean up
CLEANFILES = input.nml *.nc* *.out *_table
