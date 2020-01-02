# This is an automake file for the test_fms/diag_manager directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I${top_builddir}/mpp	\
-I${top_builddir}/fms -I${top_builddir}/diag_manager		\
-I${top_builddir}/time_manager -I${top_builddir}/constants

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_diag_manager

# This is the source code for the test.
test_diag_manager_SOURCES = test_diag_manager.F90

# Run the test.
TESTS = test_diag_manager2.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = input.nml_base diagTables test_diag_manager2.sh

CLEANFILES = input.nml *.nc *.out
