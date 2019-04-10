AM_CPPFLAGS = -I${top_builddir}/field_manager
AM_CPPFLAGS += -I${top_builddir}/mpp
check_PROGRAMS = test_field_manager
test_field_manager_SOURCES = test_field_manager.F90
test_field_manager_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_field_manager.sh
