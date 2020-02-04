# This is an automake file for the test_fms/diag_manager directory of the FMS
# package.

# J.Liptak, R.Menzel

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I${top_builddir}/mpp -I${top_builddir}/fms2_io

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_fms2_io test_atmosphere_io

# This is the source code for the test.
test_fms2_io_SOURCES = argparse.F90 test_fms2_io.F90
test_atmosphere_io_SOURCES = argparse.F90 setup.F90 test_atmosphere_io.F90
EXTRA_DIST = test_atmosphere_io.sh input_base.nml 

argparse.mod: argparse.$(OBJEXT)
setup.mod: setup.$(OBJEXT)
setup.$(OBJEXT): argparse.mod
test_atmosphere_io.$(OBJEXT): setup.mod
test_fms2_io.$(OBJEXT): argparse.mod

# Run the test program.
TESTS = test_fms2_io test_atmosphere_io.sh

# Set srcdir as evironment variable to be reference in the job script
TESTS_ENVIRONMENT = srcdir="$(top_srcdir)"

CLEANFILES = *.mod *.nc input.nml logfile.000000.out

## Build this test program.
#check_PROGRAMS = test_fms2_io
#
## This is the source code for the test.
#test_fms2_io_SOURCES = argparse.F90 test_fms2_io.F90
#test_fms2_io_LDADD = ${top_builddir}/libFMS/libFMS.la
#
#argparse.mod: argparse.$(OBJEXT)
#test_fms2_io.$(OBJEXT): argparse.mod
#
#CLEANFILES = *.mod
#
## Run the test program.
#TESTS = test_fms2_io.sh
