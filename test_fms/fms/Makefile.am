# This is an automake file for the test_fms/fms directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/fms -I${top_builddir}/mpp	\
-I${top_srcdir}/include

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_fms_io test_unstructured_fms_io

# These are the sources for the tests.
test_fms_io_SOURCES = test_fms_io.F90
test_unstructured_fms_io_SOURCES = test_unstructured_fms_io.F90

# Run the test programs.
TESTS = test_fms_io2.sh

# These will also be included in the distribution.
EXTRA_DIST = test_fms_io2.sh
