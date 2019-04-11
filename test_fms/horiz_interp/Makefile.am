AM_CPPFLAGS = -I${top_builddir}/horiz_interp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
check_PROGRAMS = test_horiz_interp
test_horiz_interp_SOURCES = test_horiz_interp.F90
test_horiz_interp_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_horiz_interp.sh
