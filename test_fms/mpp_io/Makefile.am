# This is an automake file for the test_fms/mpp directory of the FMS
# package.

# uramirez, Ed Hartnett

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_srcdir}/include -I${top_builddir}/mpp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_mpp_io

# This is the source code for the test.
test_mpp_io_SOURCES = test_mpp_io.F90

# Run the test program.
TESTS = test_mpp_io.sh

# These files will also be distributed.
EXTRA_DIST = test_mpp_io.sh test_mpp_io.bats input_base.nml

# Clean up
CLEANFILES = input.nml *.nc* *.out
