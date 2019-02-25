# This is an automake file for the drifters directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

noinst_LTLIBRARIES = libdrifters_core.la libdrifters_input.la	\
libdrifters_comm.la libdrifters_io.la libdrifters.la		\
libcloud_interpolator.la libquicksort.la

libdrifters_la_SOURCES = drifters.F90
libdrifters_core_la_SOURCES = drifters_core.F90
libdrifters_comm_la_SOURCES = drifters_comm.F90
libdrifters_input_la_SOURCES = drifters_input.F90
libdrifters_io_la_SOURCES = drifters_io.F90
libcloud_interpolator_la_SOURCES = cloud_interpolator.F90
libquicksort_la_SOURCES = quicksort.F90
