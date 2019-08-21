# This is an automake file for the test_fms/fms directory of the FMS
# package.

# Find the fms_mod.mod file.
test_fms_io_CPPFLAGS = -I${top_srcdir}/fms
test_fms_io_CPPFLAGS += -I${top_srcdir}/mpp
test_fms_io_CPPFLAGS += -I${top_srcdir}/include

test_unstructured_fms_io_CPPFLAGS = -I${top_srcdir}/fms
test_unstructured_fms_io_CPPFLAGS += -I${top_srcdir}/mpp
test_unstructured_fms_io_CPPFLAGS += -I${top_srcdir}/include

# Copy over other needed files to the srcdir
EXTRA_DIST =test_fms_io.sh
EXTRA_DIST += test_fms_io.bats
EXTRA_DIST +=  test_unstructured_fms_io.sh
EXTRA_DIST +=  test_unstructured_fms_io.bats

# Build this test program.
check_PROGRAMS = test_fms_io test_unstructured_fms_io

# This is the source code for the test.
test_fms_io_SOURCES = test_fms_io.F90
test_fms_io_LDADD = ${top_builddir}/libFMS/libFMS.la

test_unstructured_fms_io_SOURCES = test_unstructured_fms_io.F90
test_unstructured_fms_io_LDADD = ${top_builddir}/libFMS/libFMS.la

# Run the test program.
TESTS = test_fms_io.sh test_unstructured_fms_io.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"
