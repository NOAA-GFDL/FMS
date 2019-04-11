AM_CPPFLAGS = -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_interp
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/axis_utils
AM_CPPFLAGS += -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/constants

check_PROGRAMS = test_time_interp test_time_interp_external

test_time_interp_SOURCES = test_time_interp.F90
test_time_interp_LDADD = ${top_builddir}/libFMS/libFMS.la

test_time_interp_external_SOURCES = test_time_interp_external.F90
test_time_interp_external_LDADD = ${top_builddir}/libFMS/libFMS.la

TESTS = test_time_interp.sh test_time_interp_external.sh
