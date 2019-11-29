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
if ENABLE_SHELL_TESTS
TESTS = test_fft2.sh
else
test_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL)	\
$(top_srcdir)/tap-driver.sh
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"
TESTS = test_fft.test
endif

# These files are also included in the distribution.
EXTRA_DIST = test_fft2.sh test_fft.test test_fft.bats input_base.nml

# Clean up
CLEANFILES = input.nml *.out*
