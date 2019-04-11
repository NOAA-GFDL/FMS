# This is an automake file for the test_fms/mpp directory of the FMS
# package.

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

# Build this test program.
check_PROGRAMS = test_mpp test_mpp_domains test_mpp_pset

# This is the source code for the test.
test_mpp_SOURCES = test_mpp.F90
test_mpp_LDADD = ${top_builddir}/libFMS/libFMS.la

test_mpp_domains_SOURCES = test_mpp_domains.F90
test_mpp_domains_LDADD = ${top_builddir}/libFMS/libFMS.la

test_mpp_pset_SOURCES = test_mpp_pset.F90
test_mpp_pset_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_mpp.sh test_mpp_domains.sh test_mpp_pset.sh
