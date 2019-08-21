# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# Find the fms_mod.mod file.
test_fft_CPPFLAGS = -I${top_srcdir}/fft
test_fft_CPPFLAGS += -I${top_srcdir}/mpp

# Copy over other needed files to the srcdir
EXTRA_DIST = test_fft.sh
EXTRA_DIST += test_fft.bats
EXTRA_DIST += input.nml

# Build this test program.
check_PROGRAMS = test_fft

# This is the source code for the test.
test_fft_SOURCES = test_fft.F90
test_fft_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_fft.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml *.out*
