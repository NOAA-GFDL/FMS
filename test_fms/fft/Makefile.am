# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the needed mod files.
AM_CPPFLAGS = -I${top_builddir}/fft -I${top_builddir}/mpp

# Link to the FMS library.
test_fft_LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_fft

# This is the source code for the test.
test_fft_SOURCES = test_fft.F90

# Run the test program.
TESTS = test_fft2.sh

# These files are also included in the distribution.
EXTRA_DIST = test_fft2.sh input_base.nml

# Clean up
CLEANFILES = input.nml *.out*
