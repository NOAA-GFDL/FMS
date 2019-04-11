AM_CPPFLAGS = -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/time_manager
check_PROGRAMS = test_time_manager
test_time_manager_SOURCES = test_time_manager.F90
test_time_manager_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_time_manager.sh
