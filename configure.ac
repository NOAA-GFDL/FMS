#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is the main configure file for the FMS package.
# Ed Hartnett 2/21/2019
# Seth Underwood 10/2/2019

AC_PREREQ([2.69])

# Initialize with name, version, and support email address.
AC_INIT([GFDL FMS Library],
  [2022.04.00-dev],
  [gfdl.climate.model.info@noaa.gov],
  [FMS],
  [https://www.github.com/NOAA-GFDL/FMS])

# Find out about the host we're building on.
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

# Find out about the target we're building for.
AC_CANONICAL_TARGET

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])

# Set up libtool.
LT_PREREQ([2.4])
LT_INIT()

# If building on a Cray PE system, check if CRAYPE_LINK_TYPE is 'static'.  If it
# is, then disable building shared libraries.  Note, the user can still override
# this by using --enable-shared when running the configure script.
AS_IF([test x${CRAYPE_VERSION:+yes} = "xyes"],[
  AS_IF([test x${CRAYPE_LINK_TYPE} = "xstatic"],
        [AC_DISABLE_SHARED])])


# Process user optons.
AC_ARG_ENABLE([mixed-mode],
  [AS_HELP_STRING([--enable-mixed-mode],
    [Build using mixed mode.  Enables both 64-bit and 32-bit reals in Fortran. This option will be ignored if --disable-fortran-flag-setting is also given.])])
AS_IF([test ${enable_mixed_mode:-no} = no],
  [enable_mixed_mode=no],
  [enable_mixed_mode=yes])
AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi],
    [Build with MPI support.  This option will be ignored if --disable-fortran-flag-setting is also given. (Default yes)])])
AS_IF([test ${with_mpi:-yes} = yes],
  [with_mpi=yes],
  [with_mpi=no])
AC_ARG_WITH([yaml],
  [AS_HELP_STRING([--with-yaml],
    [Build with YAML support.  This option will be ignored if --disable-fortran-flag-setting is also given. (Default no)])])
AS_IF([test ${with_yaml:-no} = no],
  [with_yaml=no],
  [with_yaml=yes])
AC_ARG_ENABLE([setting-flags],
  [AS_HELP_STRING([--enable-setting-flags],
    [Allow configure to set some compiler flags.  Disabling this will also disable any other --with or --enable options that set flags, and will only use user-provided falgs. (Default yes)])])
AS_IF([test ${enable_setting_flags:-yes} = yes],
  [enable_setting_flags=yes],
  [enable_setting_flags=no])
AC_ARG_ENABLE([code-coverage],
  [AS_HELP_STRING([--enable-code-coverage],
    [Allow configure to set flags and add check-code-coverage target for code coverage analysis (Default no)])])
AS_IF([test ${enable_code_coverage:-no} = no],
  [enable_code_coverage=no],
  [enable_code_coverage=yes])
# individual mixed precision overload macros
AC_ARG_ENABLE([overload-r4],
  [AS_HELP_STRING([--enable-overload-r4],
    [Enables the OVERLOAD_R4 macro to compile with 4 byte real routine overloads. (Default no)])])
AS_IF([test ${enable_overload_r4:-no} = yes],
  [enable_overload_r4=yes],
  [enable_overload_r4=no])
AC_ARG_ENABLE([overload-c4],
  [AS_HELP_STRING([--enable-overload-c4],
    [Enables the OVERLOAD_C4 macro to compile with 4 byte complex routine overloads. (Default no)])])
AS_IF([test ${enable_overload_c4:-no} = yes],
  [enable_overload_c4=yes],
  [enable_overload_c4=no])
AC_ARG_ENABLE([overload-c8],
  [AS_HELP_STRING([--enable-overload-c8],
    [Enables the OVERLOAD_C8 macro to compile with 8 byte real routine overloads. (Default no)])])
AS_IF([test ${enable_overload_c8:-no} = yes],
  [enable_overload_c8=yes],
  [enable_overload_c8=no])
AC_ARG_ENABLE([8byte-int],
  [AS_HELP_STRING([--disable-8byte-int],
    [Enables the no_8byte_integers macro to compile with only 4 byte integer routines. (Default no)])])
AS_IF([test ${enable_8byte_int:-no} = yes],
  [enable_8byte_int=yes],
  [enable_8byte_int=no])

# user enabled testing with input files
AC_MSG_CHECKING([whether to enable tests with input files])
AC_ARG_ENABLE([test-input],
            [AS_HELP_STRING([--enable-test-input="path/to/input"],
            [Enable tests using input netcdf files, if present in the given full directory.])])
AC_MSG_RESULT([$enable_test_input])

# require path to be given
AS_IF([test "x$enable_test_input" = "xyes"],
    [AC_MSG_ERROR([Test input enabled, but no directory given with --enable-test-input=/path])],
    [AS_IF([test "x$enable_test_input" != "xno"],
        [TEST_INPUT_PATH="$enable_test_input"],[])]
)

# if set, check directory exists
AS_IF([test "$TEST_INPUT_PATH" = ""],
  [],
  [AS_IF([test -d $TEST_INPUT_PATH],[],
    [AC_MSG_ERROR([Test input enabled, but directory $TEST_INPUT_PATH not found]) ]
  )])

# substitute path for input
AC_SUBST([TEST_INPUT_PATH])

# Does the user want to build documentation?
AC_MSG_CHECKING([whether documentation should be built (requires doxygen)])
AC_ARG_ENABLE([docs],
              [AS_HELP_STRING([--enable-docs],
                              [enable building of documentation with doxygen.])])
test "x$enable_docs" = xyes || enable_docs=no
AC_MSG_RESULT([$enable_docs])
AM_CONDITIONAL(BUILD_DOCS, [test "x$enable_docs" = xyes])

# Is doxygen installed?
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN" -a "x$enable_docs" = xyes; then
   AC_MSG_ERROR([Doxygen not found but --enable-docs used.])
fi

# If building docs, process Doxyfile.in into Doxyfile.
if test "x$enable_docs" = xyes; then
   AC_CONFIG_FILES([docs/Doxyfile])
fi

# Find the C compiler.
AC_PROG_CC
AM_PROG_CC_C_O
AC_C_CONST
AC_PROG_CPP

# Find the Fortran compiler.
AC_PROG_FC
AC_FC_MODULE_FLAG
AC_FC_MODULE_EXTENSION
AC_FC_MODULE_OUTPUT_FLAG

# Find the install program.
AC_PROG_INSTALL

# Check to see if any macros must be set to enable large (>2GB) files.
AC_SYS_LARGEFILE

# C specific checks
AC_LANG_PUSH(C)
if test $with_mpi = yes; then
  # Require MPI
  # We expect users to set CC and FC to MPI compiler wrappers, or that the correct
  # CPPFLAGS, {C,FC}FLAGS, LDFLAGS and LIBS options are given.
  GX_MPI()
fi

# Require yaml
if test $with_yaml = yes; then
  AC_CHECK_HEADERS([yaml.h], [], [AC_MSG_ERROR(["Can't find the libYAML C header file.  Set CC/CPPFLAGS/CFLAGS"])])
  AC_SEARCH_LIBS([yaml_parser_initialize], [yaml], [], [AC_MSG_ERROR(["Can't find the libYAML C library.  Set CC/LDFLAGS/LIBS"])])

  #If the test pass, define use_yaml macro
  AC_DEFINE([use_yaml], [1], [This is required to use yaml parser])

  AM_CONDITIONAL([USING_YAML], true)
else
  AM_CONDITIONAL([USING_YAML], false)
fi

# Require netCDF
AC_CHECK_HEADERS([netcdf.h], [], [AC_MSG_ERROR([Can't find the netCDF C header file.  Set CPPFLAGS/CFLAGS])])
AC_SEARCH_LIBS([nc_create], [netcdf], [], [AC_MSG_ERROR([Can't find the netCDF C library.  Set LDFLAGS/LIBS])])

# Require netCDF-4 (with HDF5).
AC_MSG_CHECKING([if netCDF was built with HDF5])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#include <netcdf_meta.h>
#if !(NC_HAS_NC4)
      choke me
#endif]])], [nc_has_nc4=yes], [nc_has_nc4=no])
AC_MSG_RESULT([$nc_has_nc4])
if test $nc_has_nc4 = no; then
  AC_MSG_ERROR([NetCDF must be built with HDF5.])
fi

# Check if Linux gettid is avaiable
AC_CHECK_FUNCS([gettid], [], [])

# Check if Linux sched_getaffinity is available
AC_CHECK_FUNCS([sched_getaffinity], [], [])

# Check if the compiler needs special OpenMP flags
AC_OPENMP()
AC_LANG_POP(C)

# Fortran specific checks
AC_LANG_PUSH(Fortran)
if test $with_mpi = yes; then
  # Require MPI
  GX_MPI()
  GX_MPI_FC_LEGACY_INTERFACE()
  # Determine if a flag is required to allow external procedure argument mismatch when
  # an explicit interface does not exist
  if test ! -z "$HAVE_MPI_FC_LEGACY"; then
    GX_FC_ALLOW_ARG_MISMATCH([dnl
    FCFLAGS="$FCFLAGS $FC_ALLOW_ARG_MISMATCH_FLAG"])
  fi
fi

# check intel compiler and coverage tools are avaiable if code coverage is enabled
if test "$enable_code_coverage" = yes; then
  if test ! -z "`$FC --version | grep ifort`"; then
    AC_CHECK_PROGS([PROFMERGE],[profmerge])
    AC_CHECK_PROGS([CODECOV], [codecov])
  else
    AC_MSG_ERROR([Intel compiler and coverage tools required for coverage report])
  fi
fi
AM_CONDITIONAL(COV, [test "$enable_code_coverage" = yes])

# Require netCDF
GX_FC_CHECK_MOD([netcdf], [], [], [AC_MSG_ERROR([Can't find the netCDF Fortran module.  Set CPPFLAGS/FCFLAGS])])
GX_FORTRAN_SEARCH_LIBS([nf90_create], [netcdff], [use netcdf], [iret = nf90_create('foo.nc', 1, ncid)], [],
  [AC_MSG_ERROR([Can't find the netCDF Fortran library.  Set LDFLAGS/LIBS])])

# Check if Fortran compiler has the Class, Character array assign bug
GX_FC_CLASS_CHAR_ARRAY_BUG_CHECK()

# Check if Fortran compiler has Cray pointer support
GX_FC_CRAY_POINTER_FLAG()

# Check if Fortran compiler and system have quad precision support
GX_FC_QUAD_PRECISION()

# Check if Fortran compiler supports reading namelists from internal files
GX_FC_INTERNAL_FILE_NML()

# Check if the compiler needs special OpenMP flags
AC_OPENMP()

AC_LANG_POP(Fortran)

# We passed all the tests.  Set the required defines.
AC_DEFINE([use_netCDF], [1], [This is required for the library to build])
if test $with_mpi = yes; then
  AC_DEFINE([use_libMPI], [1], [This is required for the library to build])
fi

# Set any required compile flags.  This will not be done if the user wants to
# define all their own flags.
if test $enable_setting_flags = yes; then
  # Make sure the compiler is seeing this as free-formatted, not
  # fixed-formatted, fortran code.
  AC_FC_FREEFORM()

  # Check that long lines of Fortran code can be handled. This will add
  # necessary fortran flags.
  AC_FC_LINE_LENGTH([unlimited])

  # Will we build with default 64-bit reals in Fortran, or do mixed mode?
  if test $enable_mixed_mode = yes; then
    GX_FC_DEFAULT_REAL_KIND4_FLAG([dnl
      FCFLAGS="$FCFLAGS $FC_DEFAULT_REAL_KIND8_FLAG"])
    AC_DEFINE([OVERLOAD_R4], [1], [Set to overload the R4 Fortran routines])
    AC_DEFINE([OVERLOAD_R8], [1], [Set to overload the R8 Fortran routines])
  else
    GX_FC_DEFAULT_REAL_KIND8_FLAG([dnl
      FCFLAGS="$FCFLAGS $FC_DEFAULT_REAL_KIND8_FLAG"])
  fi
  # individual mixed precision overloads
  if test $enable_overload_r4 = yes; then
    AC_DEFINE([OVERLOAD_R4], [1], [Set to overload with the R4 Fortran routines])
  fi
  if test $enable_overload_c4 = yes; then
    AC_DEFINE([OVERLOAD_C4], [1], [Set to overload with the C4 Fortran routines])
  fi
  if test $enable_overload_c8 = yes; then
    AC_DEFINE([OVERLOAD_C8], [1], [Set to overload with the C8 Fortran routines])
  fi
  if test $enable_8byte_int = yes; then
    AC_DEFINE([no_8byte_integers], [1], [Set to disable 8 byte integer Fortran routines])
  fi

  # Add Cray Pointer support flag
  if test ! -z "$FC_CRAY_POINTER_FLAG"; then
    FCFLAGS="$FCFLAGS $FC_CRAY_POINTER_FLAG"
  fi

  # Add OpenMP flags
  if test ! -z "$OPENMP_CFLAGS"; then
    CFLAGS="$CFLAGS $OPENMP_CFLAGS"
    CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
  fi
  if test ! -z "$OPENMP_FCFLAGS"; then
    FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
  fi

  # Add code coverage flags
  if test $enable_code_coverage = yes; then
    FCFLAGS="$FCFLAGS -prof-gen=srcpos"
    CFLAGS="$CFLAGS -prof-gen=srcpos"
  fi
  # add yaml flag
  if test $with_yaml = yes; then
    LDFLAGS="$LDFLAGS -lyaml"
  fi
fi

# Find which mpi launcher to use
AC_CHECK_PROGS([MPI_LAUNCHER],[srun aprun mpirun])

# Check if the launcher can oversubscribe the MPI processes
AS_IF([$MPI_LAUNCHER --oversubscribe hostname >/dev/null 2>&1], \
      [ AC_SUBST([OVERSUBSCRIBE], [--oversubscribe])])


# Compiler with version information. This consists of the full path
# name of the compiler and the reported version number.
AC_SUBST([CC_VERSION])
# Strip anything that looks like a flag off of $CC
CC_NOFLAGS=`echo $CC | sed 's/ -.*//'`

if `echo $CC_NOFLAGS | grep ^/ >/dev/null 2>&1`; then
  CC_VERSION="$CC"
else
  CC_VERSION="$CC";
  for x in `echo $PATH | sed -e 's/:/ /g'`; do
    if test -x $x/$CC_NOFLAGS; then
      CC_VERSION="$x/$CC"
      break
    fi
  done
fi
if test -n "$cc_version_info"; then
  CC_VERSION="$CC_VERSION ( $cc_version_info)"
fi

AC_SUBST([FC_VERSION])
# Strip anything that looks like a flag off of $FC
FC_NOFLAGS=`echo $FC | sed 's/ -.*//'`

if `echo $FC_NOFLAGS | grep ^/ >/dev/null 2>&1`; then
  FC_VERSION="$FC"
else
  FC_VERSION="$FC";
  for x in `echo $PATH | sed -e 's/:/ /g'`; do
    if test -x $x/$FC_NOFLAGS; then
      FC_VERSION="$x/$FC"
      break
    fi
  done
fi
if test -n "$fc_version_info"; then
  FC_VERSION="$FC_VERSION ( $fc_version_info )"
fi

# Check if gcc is 11.1 for class(*) select type bug
AC_MSG_CHECKING([if using gcc 11.1.0])
if [ test -n "`$FC --version | grep GNU | grep 11\.1\..`" ]; then
  AC_MSG_RESULT([yes])
  AC_MSG_ERROR([Compilation with gcc and gfortran 11.1.0 is unsupported \
by this version of FMS due to a bug in the compiler. Please use a different version of gcc/gfortran.])
else
  AC_MSG_RESULT([no])
fi

#####
# Create output variables from various
# shell variables, for use in generating
# libnetcdf.settings.
#####
AC_SUBST([enable_shared])
AC_SUBST([enable_static])
AC_SUBST([CFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([NOUNDEFINED])

# Define the directory where *.mod files will reside.
# Used in Makefiles.
AC_SUBST([MODDIR],[\$\(top_builddir\)/.mods])

# These files will be created when the configure script is run.
AC_CONFIG_FILES([
  Makefile
  affinity/Makefile
  amip_interp/Makefile
  time_interp/Makefile
  time_manager/Makefile
  constants/Makefile
  constants4/Makefile
  platform/Makefile
  fms/Makefile
  fms2_io/Makefile
  mpp/Makefile
  tridiagonal/Makefile
  tracer_manager/Makefile
  topography/Makefile
  mosaic/Makefile
  mosaic2/Makefile
  monin_obukhov/Makefile
  memutils/Makefile
  interpolator/Makefile
  horiz_interp/Makefile
  field_manager/Makefile
  exchange/Makefile
  drifters/Makefile
  diag_manager/Makefile
  data_override/Makefile
  column_diagnostics/Makefile
  block_control/Makefile
  axis_utils/Makefile
  astronomy/Makefile
  coupler/Makefile
  diag_integral/Makefile
  sat_vapor_pres/Makefile
  random_numbers/Makefile
  libFMS/Makefile
  docs/Makefile
  parser/Makefile
  string_utils/Makefile
  test_fms/test-lib.sh
  test_fms/intel_coverage.sh
  test_fms/Makefile
  test_fms/diag_manager/Makefile
  test_fms/data_override/Makefile
  test_fms/exchange/Makefile
  test_fms/monin_obukhov/Makefile
  test_fms/drifters/Makefile
  test_fms/interpolator/Makefile
  test_fms/fms2_io/Makefile
  test_fms/fms/Makefile
  test_fms/mpp/Makefile
  test_fms/mpp_io/Makefile
  test_fms/time_interp/Makefile
  test_fms/time_manager/Makefile
  test_fms/horiz_interp/Makefile
  test_fms/field_manager/Makefile
  test_fms/axis_utils/Makefile
  test_fms/mosaic/Makefile
  test_fms/affinity/Makefile
  test_fms/coupler/Makefile
  test_fms/parser/Makefile
  test_fms/string_utils/Makefile
  FMS.pc
  ])

AC_OUTPUT()
