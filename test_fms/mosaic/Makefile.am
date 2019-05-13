AM_CPPFLAGS = -I${top_builddir}/mosaic
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
check_PROGRAMS = test_mosaic
test_mosaic_SOURCES = test_mosaic.F90
test_mosaic_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_mosaic.sh
