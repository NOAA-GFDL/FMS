# This is the main configure file for the FMS package.
# Ed Hartnett 2/21/2019

AC_PREREQ([2.59])

# Initialize with name, version, and support email address.
AC_INIT([FMS], [2.0-development], [])

# Find out about the host we're building on.
AC_CANONICAL_HOST

# Find out about the target we're building for.
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])

# Keep libtool macros in an m4 directory.
AC_CONFIG_MACRO_DIR([m4])

# Set up libtool.
LT_PREREQ([2.4])
LT_INIT()

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

# These defines are required for the build.
AC_DEFINE([use_netCDF], [1])
AC_DEFINE([use_libMPI], [1])

# These files will be created when the configure script is run.
AC_CONFIG_FILES([Makefile
        include/Makefile
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
        diag_integral/Makefile
        data_override/Makefile
        column_diagnostics/Makefile
        block_control/Makefile
        axis_utils/Makefile
        astronomy/Makefile
        coupler/Makefile
        sat_vapor_pres/Makefile
        random_numbers/Makefile
        libFMS/Makefile
        test_fms/Makefile
        ])
AC_OUTPUT()
