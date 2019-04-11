AM_CPPFLAGS = -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/diag_manager
AM_CPPFLAGS += -I${top_builddir}/interpolator
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/time_interp
check_PROGRAMS = test_interpolator
test_interpolator_SOURCES = test_interpolator.F90
test_interpolator_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_interpolator.sh

