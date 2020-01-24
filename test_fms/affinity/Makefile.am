# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez

# Find the fms and mpp mod files.
AM_CPPFLAGS = -I${top_builddir}/mpp -I${top_builddir}/affinity

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_affinity

# This is the source code for the test.
test_affinity_SOURCES = test_affinity.F90

# Run the test program.
TESTS = test_affinity2.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = input_base.nml test_affinity2.sh

# Clean up
CLEANFILES = input.nml *.out*
