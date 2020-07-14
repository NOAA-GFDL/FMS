# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I${top_builddir}/mpp	\
-I${top_builddir}/fms -I${top_builddir}/constants		\
-I${top_builddir}/time_manager -I${top_builddir}/diag_manager	\
-I${top_builddir}/data_override

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_data_override

# This is the source code for the test.
test_data_override_SOURCES = test_data_override.F90

# Run the test program.
TESTS = test_data_override2.sh

# Include these files with the distribution.
EXTRA_DIST = input_base.nml diag_table_base data_table_base	\
test_data_override2.sh

# Clean up
CLEANFILES = input.nml *.nc* *.out diag_table data_table
