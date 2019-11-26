# This is the main configure file for the FMS package.
# Ed Hartnett 2/21/2019

AC_PREREQ([2.59])

# Initialize with name, version, and support email address.
AC_INIT([GFDL FMS Library],
  [2019.01],
  [gfdl.climate.model.info@noaa.gov],
  [FMS],
  [https://www.gfdl.noaa.gov/fms])

# Find out about the host we're building on.
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

# Find out about the target we're building for.
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])
AC_REQUIRE_AUX_FILE([tap-driver.sh])

# Keep libtool macros in an m4 directory.
AC_CONFIG_MACRO_DIR([m4])

# Set up libtool.
LT_PREREQ([2.4])
LT_INIT()

# Process user optons.

# Does the user want to turn off all attempts to set fortran flags?
# This is necessary to handle situations where the configure-provided
# flags will be wrong for a particular install. With this setting, the
# user is in complete control of Fortran flags.
AC_MSG_CHECKING([whether we should attempt to set Fortran flags for the user])
AC_ARG_ENABLE([fortran-flag-setting],
  [AS_HELP_STRING([--disable-fortran-flag-setting],
    [Disable any attempts by configure to set fortran flags. With this option, only user-provided flags (set in FCFLAGS) will be used compiling Fortran code.])])
test "x$enable_fortran_flag_setting" = xno || enable_fortran_flag_setting=yes
AC_MSG_RESULT([$enable_fortran_flag_setting])

# Find the C compiler.
AC_PROG_CC
AM_PROG_CC_C_O
AC_C_CONST
AC_PROG_CPP

# Find the Fortran compiler.
AC_PROG_FC
AC_PROG_F77

# Find the install program.
AC_PROG_INSTALL

# Is bats installed? It is required for testing.
AC_CHECK_PROGS([FMS_BATS], [bats])
if test -z "$FMS_BATS"; then
   AC_MSG_ERROR([Cannot find bats utility. Install bats and try again.])
fi

# Check to see if any macros must be set to enable large (>2GB) files.
AC_SYS_LARGEFILE

# Require MPI.
AC_CHECK_FUNC([MPI_Init], [], [AC_MSG_ERROR([MPI C library required to build FMS])])

# Check for netCDF C library.
AC_SEARCH_LIBS([nc_create], [netcdf], [],
                            [AC_MSG_ERROR([Can't find or link to the netcdf C library, set CPPFLAGS/LDFLAGS.])])

# Check for netCDF Fortran library.
AC_LANG_PUSH(Fortran)
AC_SEARCH_LIBS([nf_create], [netcdff], [],
                            [AC_MSG_ERROR([Can't find or link to the netcdf Fortran library, set CPPFLAGS/LDFLAGS.])])
AC_LANG_POP(Fortran)

# Require netCDF.
#AC_CHECK_FUNC([nf_open], [], [AC_MSG_ERROR([NetCDF Fortran library required to build FMS])])

# Do we want to try and set Fortran flags for the user?
if test $enable_fortran_flag_setting = "yes"; then
   # Make sure the compiler is seeing this as free-formatted, not
   # fixed-formatted, fortran code.
   AC_FC_FREEFORM()

   # Check that long lines of Fortran code can be handled. This will add
   # necessary fortran flags.
   AC_FC_LINE_LENGTH([unlimited])
fi

# These defines are required for the build.
AC_DEFINE([use_netCDF], [1])
AC_DEFINE([use_libMPI], [1])

# Define __APPLE__ macro if on a Apple/Mac OS X
AS_IF([test "x$target_vendor" = "xapple"], [AC_DEFINE([__APPLE__], [ ])])

# Define an AM_CONDITIONAL to determine if you are on a CRAY
AM_CONDITIONAL([CRAY], [test `env | grep CRAY | wc -l` -gt 0])

# Define an AM_CONDITIONAL to determine if you are on an Apple/Mac OS X
AM_CONDITIONAL([APPLE], [test "x$target_vendor" = "xapple"])

# These files will be created when the configure script is run.
AC_CONFIG_FILES([Makefile
        include/Makefile
        affinity/Makefile
        amip_interp/Makefile
        time_interp/Makefile
        time_manager/Makefile
        constants/Makefile
        platform/Makefile
        fms/Makefile
        mpp/Makefile
        mpp/include/Makefile
        tridiagonal/Makefile
        tracer_manager/Makefile
        topography/Makefile
        station_data/Makefile
        oda_tools/Makefile
        mosaic/Makefile
        monin_obukhov/Makefile
        memutils/Makefile
        interpolator/Makefile
        horiz_interp/Makefile
        field_manager/Makefile
        fft/Makefile
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
        test_fms/Makefile
        test_fms/diag_manager/Makefile
        test_fms/data_override/Makefile
        test_fms/exchange/Makefile
        test_fms/monin_obukhov/Makefile
        test_fms/drifters/Makefile
        test_fms/interpolator/Makefile
        test_fms/fms/Makefile
        test_fms/mpp/Makefile
        test_fms/mpp_io/Makefile
        test_fms/time_interp/Makefile
        test_fms/time_manager/Makefile
        test_fms/horiz_interp/Makefile
        test_fms/field_manager/Makefile
        test_fms/fft/Makefile
        test_fms/axis_utils/Makefile
        test_fms/mosaic/Makefile
        ])
AC_OUTPUT()
