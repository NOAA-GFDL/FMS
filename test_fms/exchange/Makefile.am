# This is an automake file for the test_fms/exchange directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_xgrid_CPPFLAGS = -I${top_srcdir}/include
test_xgrid_CPPFLAGS += -I${top_builddir}/mpp
test_xgrid_CPPFLAGS += -I${top_builddir}/fms
test_xgrid_CPPFLAGS += -I${top_builddir}/constants
test_xgrid_CPPFLAGS += -I${top_builddir}/exchange
test_xgrid_CPPFLAGS += -I${top_builddir}/mosaic
test_xgrid_CPPFLAGS += -I${top_builddir}/coupler

# Copy over other needed files to the srcdir
EXTRA_DIST = test_xgrid.sh
EXTRA_DIST += input.nml
#EXTRA_DIST += INPUT
EXTRA_DIST += test_xgrid.bats

# Build this test program.
check_PROGRAMS = test_xgrid

# This is the source code for the test.
test_xgrid_SOURCES = test_xgrid.F90
test_xgrid_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_xgrid.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"
TESTS_ENVIRONMENT += builddir="$(top_builddir)"

# Clean up
DISTCLEANFILES = input.nml
DISTCLEANFILES += *.nc
DISTCLEANFILES += *.out
