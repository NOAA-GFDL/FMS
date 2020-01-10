# This is an automake file for the test_fms/mpp_io directory of the
# FMS package.

# uramirez, Ed Hartnett

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_srcdir}/include -I${top_builddir}/mpp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_mpp_io

# This is the source code for the test.
test_mpp_io_SOURCES = test_mpp_io.F90

# Run the test program.
TESTS = test_mpp_io2.sh

# These files will also be distributed.
EXTRA_DIST = test_mpp_io2.sh input_base.nml

# Clean up
CLEANFILES = input.nml *.nc* *.out
