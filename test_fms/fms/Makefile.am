AM_CPPFLAGS = -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/include

check_PROGRAMS = test_fms_io test_unstructured_fms_io

test_fms_io_SOURCES = test_fms_io.F90
test_fms_io_LDADD = ${top_builddir}/libFMS/libFMS.la

test_unstructured_fms_io_SOURCES = test_unstructured_fms_io.F90
test_unstructured_fms_io_LDADD = ${top_builddir}/libFMS/libFMS.la

TESTS = test_fms_io.sh test_unstructured_fms_io.sh
