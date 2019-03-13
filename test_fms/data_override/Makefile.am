# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/data_override

# Build this test program.
check_PROGRAMS = test_data_override

# This is the source code for the test.
test_data_override_SOURCES = test_data_override.F90
test_data_override_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_data_override.sh
