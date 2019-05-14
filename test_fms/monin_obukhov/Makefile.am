AM_CPPFLAGS = -I${top_builddir}/monin_obukhov
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/mpp
check_PROGRAMS = test_monin_obukhov
test_monin_obukhov_SOURCES = test_monin_obukhov.F90
test_monin_obukhov_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_monin_obukhov.sh 
