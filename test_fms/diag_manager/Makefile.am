# This is an automake file for the test_fms/diag_manager directory of the FMS
# package.

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/include

# Build this test program.
check_PROGRAMS = test_diag_manager

# This is the source code for the test.
test_diag_manager_SOURCES = test_diag_manager.F90
test_diag_manager_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_diag_manager.sh
