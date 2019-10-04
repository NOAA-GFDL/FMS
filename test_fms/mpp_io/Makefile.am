# This is an automake file for the test_fms/mpp directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_mpp_io_CPPFLAGS = -I${top_srcdir}/include
test_mpp_io_CPPFLAGS += -I${top_builddir}/mpp

# Build this test program.
check_PROGRAMS = test_mpp_io

# This is the source code for the test.
test_mpp_io_SOURCES = test_mpp_io.F90
test_mpp_io_LDADD = ${top_builddir}/libFMS/libFMS.la

# Copy over other needed files to the srcdir
EXTRA_DIST = test_mpp_io.sh
EXTRA_DIST += test_mpp_io.bats
EXTRA_DIST += input.nml

# Run the test program.
TESTS = test_mpp_io.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = ${top_builddir}/test_fms/mpp_io/input.nml
CLEANFILES += ${top_builddir}/test_fms/mpp_io/*.nc*
CLEANFILES += ${top_builddir}/test_fms/mpp_io/*.out
