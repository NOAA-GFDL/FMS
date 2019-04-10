AM_CPPFLAGS = -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
check_PROGRAMS = test_axis_utils
test_axis_utils_SOURCES = test_axis_utils.F90
test_axis_utils_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_axis_utils.sh
