# This is an automake file for the test_fms directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
LDADD = ${top_builddir}/libFMS/.libs/libFMS.la

# Build this test program.
check_PROGRAMS = tst_fms1

# This is the source code for the test.
tst_fms1_SOURCES = tst_fms1.F90

# Run the test program.
TESTS = tst_fms1
