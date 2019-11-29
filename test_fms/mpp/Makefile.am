# This is an automake file for the test_fms/mpp directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the needed mod and inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I${top_builddir}/mpp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_mpp test_mpp_domains test_mpp_pset

# These are the sources for the tests.
test_mpp_SOURCES = test_mpp.F90
test_mpp_domains_SOURCES = test_mpp_domains.F90
test_mpp_pset_SOURCES = test_mpp_pset.F90

# Run the test programs.
if ENABLE_SHELL_TESTS
TESTS = test_mpp_domains2.sh test_mpp_pset2.sh test_mpp2.sh
else
test_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) $(top_srcdir)/tap-driver.sh
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"
TESTS = test_mpp_domains.test test_mpp_pset.test test_mpp.test
endif

# These files will also be included in the distribution.
EXTRA_DIST = input_base.nml test_mpp.test test_mpp.bats \
test_mpp_domains.test test_mpp_domains.bats test_mpp_pset.test \
test_mpp_pset.bats test_mpp_domains2.sh test_mpp_pset2.sh test_mpp2.sh

# Clean up
CLEANFILES = input.nml *.out*
