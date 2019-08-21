# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# Find the fms_mod.mod file.
test_field_manager_CPPFLAGS = -I${top_srcdir}/field_manager
test_field_manager_CPPFLAGS += -I${top_srcdir}/mpp

# Copy over other needed files to the srcdir
EXTRA_DIST = test_field_manager.sh
EXTRA_DIST +=  test_field_manager.bats
EXTRA_DIST += input.nml
EXTRA_DIST += field_table

# Build this test program.
check_PROGRAMS = test_field_manager

# This is the source code for the test.
test_field_manager_SOURCES = test_field_manager.F90
test_field_manager_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_field_manager.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml *.out* field_table
