# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# Find the fms_mod.mod file.
test_horiz_interp_CPPFLAGS = -I${top_srcdir}/horiz_interp
test_horiz_interp_CPPFLAGS += -I${top_srcdir}/fms
test_horiz_interp_CPPFLAGS += -I${top_srcdir}/constants
test_horiz_interp_CPPFLAGS += -I${top_srcdir}/mpp

test2_horiz_interp_CPPFLAGS = -I${top_srcdir}/horiz_interp
test2_horiz_interp_CPPFLAGS += -I${top_srcdir}/fms
test2_horiz_interp_CPPFLAGS += -I${top_srcdir}/constants
test2_horiz_interp_CPPFLAGS += -I${top_srcdir}/mpp
test2_horiz_interp_CPPFLAGS += -I${top_srcdir}/mosaic
test2_horiz_interp_CPPFLAGS += -I${top_srcdir}/axis_utils

# Copy over other needed files to the srcdir
EXTRA_DIST = test_horiz_interp.sh
EXTRA_DIST +=  test_horiz_interp.bats
EXTRA_DIST += test2_horiz_interp.sh
EXTRA_DIST +=  test2_horiz_interp.bats
EXTRA_DIST += input.nml

# Build this test program.
check_PROGRAMS = test_horiz_interp test2_horiz_interp

# This is the source code for the test.
test_horiz_interp_SOURCES = test_horiz_interp.F90
test_horiz_interp_LDADD = ${top_builddir}/libFMS/libFMS.la

test2_horiz_interp_SOURCES = test2_horiz_interp.F90
test2_horiz_interp_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_horiz_interp.sh test2_horiz_interp.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml *.out* 
