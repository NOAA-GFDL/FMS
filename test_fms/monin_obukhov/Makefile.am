# This is an automake file for the test_fms/monin_obukhov directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_monin_obukhov_CPPFLAGS = -I${top_builddir}/monin_obukhov
test_monin_obukhov_CPPFLAGS += -I${top_builddir}/fms
test_monin_obukhov_CPPFLAGS += -I${top_builddir}/mpp

# Copy over other needed files to the srcdir
EXTRA_DIST = input_base.nml
EXTRA_DIST += test_monin_obukhov.sh
EXTRA_DIST += test_monin_obukhov.bats

# Build this test program.
check_PROGRAMS = test_monin_obukhov

# This is the source code for the test.
test_monin_obukhov_SOURCES = test_monin_obukhov.F90
test_monin_obukhov_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_monin_obukhov.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml *.out
