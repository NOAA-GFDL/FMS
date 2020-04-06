# This is an automake file for the test_fms/fms2_io directory of the FMS
# package.

# J.Liptak, R.Menzel

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/mpp -I${top_builddir}/fms2_io

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_fms2_io test_atmosphere_io

# This is the source code for the test.
test_fms2_io_SOURCES = argparse.F90 test_fms2_io.F90 create_atmosphere_domain.inc \
                       create_land_domain.inc create_ocean_domain.inc \
                       ocean_restart_file_test.inc land_compressed_restart_file_test.inc \
                       land_unstructured_restart_file_test.inc
test_atmosphere_io_SOURCES = argparse.F90 setup.F90 test_atmosphere_io.F90 \
                             atmosphere_restart_file_test.inc

EXTRA_DIST = test_fms2_io.sh test_atmosphere_io.sh

argparse.mod: argparse.$(OBJEXT)
setup.mod: setup.$(OBJEXT)
setup.$(OBJEXT): argparse.mod
test_atmosphere_io.$(OBJEXT): setup.mod
test_fms2_io.$(OBJEXT): argparse.mod

# Run the test program.
TESTS = test_fms2_io.sh test_atmosphere_io.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

CLEANFILES = *.mod *.nc *.nc.* input.nml logfile.000000.out
