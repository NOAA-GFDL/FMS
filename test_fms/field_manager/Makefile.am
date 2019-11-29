# This is an automake file for the test_fms/data_override directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the needed mod files.
AM_CPPFLAGS = -I${top_builddir}/field_manager -I${top_builddir}/mpp

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_field_manager

# This is the source code for the test.
test_field_manager_SOURCES = test_field_manager.F90

# Run the test program.
if ENABLE_SHELL_TESTS
TESTS = test_field_manager2.sh
else
TESTS = test_field_manager.test
test_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL)	\
$(top_srcdir)/tap-driver.sh
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"
endif

# These files will also be included in the distribution.
EXTRA_DIST = test_field_manager2.sh test_field_manager.test	\
test_field_manager.bats input_base.nml field_table_base

# Clean up
CLEANFILES = input.nml *.out* field_table
