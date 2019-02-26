# This is an automake file for the libFMS directory of the FMS
# package. This is the final packaging of the library.

# Ed Hartnett 2/22/19

lib_LTLIBRARIES = libFMS.la

# These linker flags specify libtool version info.
# See http://www.gnu.org/software/libtool/manual/libtool.html#Libtool-versioning
# for information regarding incrementing `-version-info`.
libFMS_la_LDFLAGS = -version-info 2:0:0

libFMS_la_LIBADD = ${top_builddir}/constants/libconstants.la
libFMS_la_LIBADD += ${top_builddir}/topography/libtopography.la
libFMS_la_LIBADD += ${top_builddir}/coupler/libatmos_ocean_fluxes.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libnsclock.la
libFMS_la_LIBADD += ${top_builddir}/oda_tools/liboda_types.la
libFMS_la_LIBADD += ${top_builddir}/time_manager/libget_cal_time.la
libFMS_la_LIBADD += ${top_builddir}/monin_obukhov_kernel/libmonin_obukhov_kernel.la
libFMS_la_LIBADD += ${top_builddir}/fft/libfft.la
libFMS_la_LIBADD += ${top_builddir}/amip_interp/libamip_interp.la

# At least one source file must be included.
libFMS_la_SOURCES = ${top_builddir}/include/file_version.h
