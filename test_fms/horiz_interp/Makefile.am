# This is an automake file for the test_fms/horiz_interp directory of
# the FMS package.

# uramirez, Ed Hartnett

# Find the needed mod files.
AM_CPPFLAGS = -I${top_builddir}/horiz_interp -I${top_builddir}/fms	\
-I${top_builddir}/constants -I${top_builddir}/mpp			\
-I${top_builddir}/mosaic -I${top_builddir}/axis_utils

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_horiz_interp test2_horiz_interp

# These are the sources for the tests.
test_horiz_interp_SOURCES = test_horiz_interp.F90
test2_horiz_interp_SOURCES = test2_horiz_interp.F90

# Run the test programs.
TESTS = test_horiz_interp2.sh

# These files will also be included in the distribution.
EXTRA_DIST = input_base.nml test_horiz_interp2.sh

# Clean up
CLEANFILES = input.nml *.out*
