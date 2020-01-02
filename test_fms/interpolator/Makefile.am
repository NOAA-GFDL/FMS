# This is an automake file for the test_fms/interpolator directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/mpp -I${top_builddir}/fms	\
-I${top_builddir}/time_manager -I${top_builddir}/diag_manager	\
-I${top_builddir}/interpolator -I${top_builddir}/constants	\
-I${top_builddir}/time_interp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_interpolator

# This is the source code for the test.
test_interpolator_SOURCES = test_interpolator.F90

# Run the test program.
TESTS = test_interpolator2.sh

# These files are also included in the distribution.
EXTRA_DIST = test_interpolator2.sh input_base.nml diag_table_base

# Clean up
CLEANFILES = input.nml diag_table
