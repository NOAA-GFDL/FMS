
AM_CPPFLAGS = -I${top_builddir}/drifters
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/include
AM_CPPFLAGS += -I${top_builddir}/fms

check_PROGRAMS = test_cloud_interpolator test_drifters_input test_drifters_comm test_drifters_core test_drifters_io test_quicksort test_drifters

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

TESTS = test_drifters_io.sh test_cloud_interpolator.sh test_drifters_input.sh test_drifters_comm.sh test_drifters_core.sh test_quicksort.sh test_drifters.sh
