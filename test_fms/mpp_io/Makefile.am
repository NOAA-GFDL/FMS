# This is an automake file for the test_fms/mpp directory of the FMS
# package.

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build this test program.
check_PROGRAMS = test_mpp_io

# This is the source code for the test.
test_mpp_io_SOURCES = test_mpp_io.F90
test_mpp_io_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_mpp_io.sh
