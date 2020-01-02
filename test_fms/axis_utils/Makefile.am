# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the fms and mpp mod files.
AM_CPPFLAGS = -I${top_builddir}/axis_utils -I${top_builddir}/fms	\
-I${top_builddir}/mpp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_axis_utils

# This is the source code for the test.
test_axis_utils_SOURCES = test_axis_utils.F90

# Run the test program.
TESTS = test_axis_utils2.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = input_base.nml test_axis_utils2.sh

# Clean up
CLEANFILES = input.nml *.out*
