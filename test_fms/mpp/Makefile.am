# This is an automake file for the test_fms/mpp directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_mpp_CPPFLAGS = -I${top_srcdir}/include
test_mpp_CPPFLAGS += -I${top_builddir}/mpp

test_mpp_domains_CPPFLAGS = -I${top_srcdir}/include
test_mpp_domains_CPPFLAGS += -I${top_builddir}/mpp

test_mpp_pset_CPPFLAGS = -I${top_srcdir}/include
test_mpp_pset_CPPFLAGS += -I${top_builddir}/mpp

# Build this test program.
check_PROGRAMS = test_mpp test_mpp_domains test_mpp_pset

# Copy over other needed files to the srcdir
EXTRA_DIST = input.nml_base
EXTRA_DIST += test_mpp.sh
EXTRA_DIST +=  test_mpp.bats
EXTRA_DIST += test_mpp_domains.sh
EXTRA_DIST +=  test_mpp_domains.bats
EXTRA_DIST += test_mpp_pset.sh
EXTRA_DIST +=  test_mpp_pset.bats

# This is the source code for the test.
test_mpp_SOURCES = test_mpp.F90
test_mpp_LDADD = ${top_builddir}/libFMS/libFMS.la

test_mpp_domains_SOURCES = test_mpp_domains.F90
test_mpp_domains_LDADD = ${top_builddir}/libFMS/libFMS.la

test_mpp_pset_SOURCES = test_mpp_pset.F90
test_mpp_pset_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_mpp_domains.sh test_mpp_pset.sh test_mpp.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = input.nml* *.out*
