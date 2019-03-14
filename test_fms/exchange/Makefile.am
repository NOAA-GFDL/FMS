# This is an automake file for the test_fms/exchange directory of the FMS
# package.

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/exchange
AM_CPPFLAGS += -I${top_builddir}/mosaic
AM_CPPFLAGS += -I${top_builddir}/coupler

# Build this test program.
check_PROGRAMS = test_xgrid

# This is the source code for the test.
test_xgrid_SOURCES = test_xgrid.F90
test_xgrid_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_xgrid.sh
