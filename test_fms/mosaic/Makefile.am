# This is an automake file for the test_fms/mosaic directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_mosaic_CPPFLAGS = -I${top_builddir}/mosaic
test_mosaic_CPPFLAGS += -I${top_builddir}/fms
test_mosaic_CPPFLAGS += -I${top_srcdir}/include
test_mosaic_CPPFLAGS += -I${top_builddir}/mpp

# Build this test program.
check_PROGRAMS = test_mosaic

# Copy over other needed files to the srcdir
EXTRA_DIST = test_mosaic.sh
EXTRA_DIST +=  test_mosaic.bats
EXTRA_DIST += input.nml

# This is the source code for the test.
test_mosaic_SOURCES = test_mosaic.F90
test_mosaic_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_mosaic.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"
