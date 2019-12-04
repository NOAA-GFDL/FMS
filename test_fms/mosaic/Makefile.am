# This is an automake file for the test_fms/mosaic directory of the
# FMS package.

# uramirez, Ed Hartnett

# Find the needed mod and include files.
AM_CPPFLAGS = -I${top_builddir}/mosaic -I${top_builddir}/fms	\
-I${top_srcdir}/include -I${top_builddir}/mpp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_mosaic

# This is the source code for the test.
test_mosaic_SOURCES = test_mosaic.F90

# Run the test program.
TESTS = test_mosaic2.sh

# These files are also included in the distribution.
EXTRA_DIST = test_mosaic2.sh input.nml
