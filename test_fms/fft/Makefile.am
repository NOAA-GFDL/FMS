AM_CPPFLAGS = -I${top_builddir}/fft
AM_CPPFLAGS += -I${top_builddir}/mpp
check_PROGRAMS = test_fft
test_fft_SOURCES = test_fft.F90
test_fft_LDADD = ${top_builddir}/libFMS/libFMS.la
TESTS = test_fft.sh
