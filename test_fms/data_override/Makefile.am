# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# Find the fms_mod.mod file.
test_data_override_CPPFLAGS = -I${top_srcdir}/include
test_data_override_CPPFLAGS += -I${top_srcdir}/mpp
test_data_override_CPPFLAGS += -I${top_srcdir}/fms
test_data_override_CPPFLAGS += -I${top_srcdir}/constants
test_data_override_CPPFLAGS += -I${top_srcdir}/time_manager
test_data_override_CPPFLAGS += -I${top_srcdir}/diag_manager
test_data_override_CPPFLAGS += -I${top_srcdir}/data_override

# Copy over other needed files to the srcdir
EXTRA_DIST = test_data_override.sh
EXTRA_DIST +=  test_data_override.bats
EXTRA_DIST += input.nml_base
EXTRA_DIST += diag_table
EXTRA_DIST += data_table
#EXTRA_DIST += INPUT

# Build this test program.
check_PROGRAMS = test_data_override

# This is the source code for the test.
test_data_override_SOURCES = test_data_override.F90
test_data_override_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_data_override.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"
TESTS_ENVIRONMENT += builddir="$(top_builddir)"

# Clean up
DISTCLEANFILES = input.nml
DISTCLEANFILES += *.nc*
DISTCLEANFILES += *.out
DISTCLEANFILES += *_table
