# This is an automake file for the test_fms/drifters directory of the FMS
# package.

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

# Find the fms_mod.mod file.
test_cloud_interpolator_CPPFLAGS = -I${top_builddir}/drifters
test_cloud_interpolator_CPPFLAGS += -I${top_builddir}/mpp
test_cloud_interpolator_CPPFLAGS += -I${top_srcdir}/include

test_drifters_io_CPPFLAGS = -I${top_builddir}/drifters
test_drifters_io_CPPFLAGS += -I${top_builddir}/mpp
test_drifters_io_CPPFLAGS += -I${top_builddir}/include

test_drifters_input_CPPFLAGS = -I${top_builddir}/drifters
test_drifters_input_CPPFLAGS += -I${top_builddir}/mpp
test_drifters_input_CPPFLAGS += -I${top_builddir}/fms
test_drifters_input_CPPFLAGS += -I${top_builddir}/include

test_drifters_comm_CPPFLAGS = -I${top_builddir}/drifters
test_drifters_comm_CPPFLAGS += -I${top_builddir}/mpp
test_drifters_comm_CPPFLAGS += -I${top_builddir}/include

test_drifters_core_CPPFLAGS = -I${top_builddir}/drifters
test_drifters_core_CPPFLAGS += -I${top_builddir}/mpp
test_drifters_core_CPPFLAGS += -I${top_builddir}/include
test_drifters_core_CPPFLAGS += -I${top_builddir}/fms

test_quicksort_CPPFLAGS = -I${top_builddir}/drifters

test_drifters_CPPFLAGS = -I${top_builddir}/drifters
test_drifters_CPPFLAGS += -I${top_builddir}/mpp
test_drifters_CPPFLAGS += -I${top_builddir}/include

# Copy over other needed files to the srcdir
EXTRA_DIST = input_base.nml
EXTRA_DIST += test_cloud_interpolator.bats
EXTRA_DIST += test_cloud_interpolator.sh
EXTRA_DIST += test_drifters_io.bats
EXTRA_DIST += test_drifters_io.sh
EXTRA_DIST += test_drifters_input.bats
EXTRA_DIST += test_drifters_input.sh
EXTRA_DIST += test_drifters_comm.bats
EXTRA_DIST += test_drifters_comm.sh
EXTRA_DIST += test_drifters_core.bats
EXTRA_DIST += test_drifters_core.sh
EXTRA_DIST += test_quicksort.bats
EXTRA_DIST += test_quicksort.sh
EXTRA_DIST += test_drifters.bats
EXTRA_DIST += test_drifters.sh
EXTRA_DIST += drifters_inp_test_3d.cdl

# Build this test program.
check_PROGRAMS = test_cloud_interpolator test_drifters_input test_drifters_comm test_drifters_core test_drifters_io test_quicksort test_drifters

# This is the source code for the tests.
test_cloud_interpolator_SOURCES = test_cloud_interpolator.F90
test_cloud_interpolator_LDADD = ${top_builddir}/libFMS/libFMS.la

test_drifters_io_SOURCES = test_drifters_io.F90
test_drifters_io_LDADD = ${top_builddir}/libFMS/libFMS.la

test_drifters_input_SOURCES = test_drifters_input.F90
test_drifters_input_LDADD = ${top_builddir}/libFMS/libFMS.la

test_drifters_comm_SOURCES = test_drifters_comm.F90
test_drifters_comm_LDADD = ${top_builddir}/libFMS/libFMS.la

test_drifters_core_SOURCES = test_drifters_core.F90
test_drifters_core_LDADD = ${top_builddir}/libFMS/libFMS.la

test_quicksort_SOURCES = test_quicksort.F90
test_quicksort_LDADD = ${top_builddir}/libFMS/libFMS.la

test_drifters_SOURCES = test_drifters.F90
test_drifters_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_drifters_io.sh test_cloud_interpolator.sh test_drifters_input.sh test_drifters_comm.sh test_drifters_core.sh test_quicksort.sh test_drifters.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

# Clean up
CLEANFILES = *out *__genmod.mod *__genmod.f90 input.nml *.nc
