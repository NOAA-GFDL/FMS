# This is an automake file for the test_fms/time_manager directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/mpp -I${top_builddir}/fms	\
-I${top_builddir}/constants -I${top_builddir}/time_manager

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_time_manager

# This is the source code for the test.
test_time_manager_SOURCES = test_time_manager.F90

# Run the test program.
TESTS = test_time_manager2.sh

# These files will be included in the distribution.
EXTRA_DIST = test_time_manager2.sh input_base.nml

# Clean up
CLEANFILES = input.nml *.out*
