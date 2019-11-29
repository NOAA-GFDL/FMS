# This is an automake file for the test_fms/drifters directory of the FMS
# package.

# uramirez, Ed Hartnett

test_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL)	\
$(top_srcdir)/tap-driver.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_builddir}/drifters -I${top_builddir}/mpp \
-I${top_srcdir}/include -I${top_builddir}/fms

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_cloud_interpolator test_drifters_input		\
test_drifters_comm test_drifters_core test_drifters_io test_quicksort	\
test_drifters

# This is the source code for the tests.
test_cloud_interpolator_SOURCES = test_cloud_interpolator.F90
test_drifters_io_SOURCES = test_drifters_io.F90
test_drifters_input_SOURCES = test_drifters_input.F90
test_drifters_comm_SOURCES = test_drifters_comm.F90
test_drifters_core_SOURCES = test_drifters_core.F90
test_quicksort_SOURCES = test_quicksort.F90
test_drifters_SOURCES = test_drifters.F90

# Run the test programs.
if ENABLE_SHELL_TESTS
TESTS = test_drifters2.sh
else
TESTS = test_drifters_io.test test_cloud_interpolator.test		\
test_drifters_comm.test test_drifters_core.test test_quicksort.test
#test_drifters.test test_drifters_input.test
endif

# Include these files in the distribution.
EXTRA_DIST = input_base.nml test_cloud_interpolator.bats	\
test_cloud_interpolator.test test_drifters_io.bats		\
test_drifters_io.test test_drifters_input.bats			\
test_drifters_input.test test_drifters_comm.bats		\
test_drifters_comm.test test_drifters_core.bats			\
test_drifters_core.test test_quicksort.bats test_quicksort.test	\
test_drifters.bats test_drifters.test drifters_inp_test_3d.cdl	\
test_drifters2.sh

# Clean up
CLEANFILES = *out *__genmod.mod *__genmod.f90 input.nml *.nc
