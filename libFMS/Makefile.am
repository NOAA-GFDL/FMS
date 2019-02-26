# This is an automake file for the libFMS directory of the FMS
# package. This is the final packaging of the library.

# Ed Hartnett 2/22/19

lib_LTLIBRARIES = libFMS.la

# These linker flags specify libtool version info.
# See http://www.gnu.org/software/libtool/manual/libtool.html#Libtool-versioning
# for information regarding incrementing `-version-info`.
libFMS_la_LDFLAGS = -version-info 2:0:0

libFMS_la_LIBADD = ${top_builddir}/constants/libconstants.la

# At least one source file must be included.
libFMS_la_SOURCES = ${top_builddir}/include/file_version.h
