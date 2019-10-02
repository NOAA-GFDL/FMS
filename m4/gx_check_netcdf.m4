# GX_CHECK_NETCDF verifies the netCDF header files (netcdf.h for C/C++ and
# netcdf.inc for Fortran 77/Fortran and Fortran module netcdf) can be found
# and used.  Also checks if the netCDF libraries (-lnetcdf for C/C++ and
# -lnetcdff for Fortran 77/Fortran) are usable.  If the netCDF header files
# are found, the function will define HAVE_header_file (in all caps).  If
# the libraries are found and usable, LIBS will be set to include the libraries.
AC_DEFUN([GX_CHECK_NETCDF], [
  AC_LANG_CASE(
    [C],[
      AC_REQUIRE([AC_PROG_CC])
      AC_CHECK_HEADERS([netcdf.h])
      AC_SEARCH_LIBS([nc_create], [netcdf])
    ],
    [C++],[
      AC_REQUIRE([AC_PROG_CXX])
      AC_CHECK_HEADERS([netcdf.h])
      AC_SEARCH_LIBS([nc_create], [netcdf])
    ],
    [Fortran 77],[
      AC_REQUIRE([AC_PROG_F77])
      dnl The include information may be included in the CPPFLAGS.  Autoconf
      dnl does not pass this flag to the Fortran compiler during the compile
      dnl or link tests.  Doing so manually.
      _gx_check_netcdf_save_FFLAGS=$FFLAGS
      FFLAGS="$CPPFLAGS $FFLAGS"
      _GX_CHECK_NETCDF_INC()
      _GX_SEARCH_NETCDFF_LIBS()
      dnl Restore FFLAGS
      FFLAGS="$_gx_check_netcdf_save_FFLAGS"
      unset _gx_check_netcdf_save_FFLAGS
    ],
    [Fortran],[
      AC_REQUIRE([AC_PROG_FC])
      dnl The include information may be included in the CPPFLAGS.  Autoconf
      dnl does not pass this flag to the Fortran compiler during the compile
      dnl or link tests.  Doing so manually.
      _gx_check_netcdf_save_FCFLAGS=$FCFLAGS
      FCFLAGS="$CPPFLAGS $FCFLAGS"
      _GX_CHECK_NETCDF_INC()
      _GX_SEARCH_NETCDFF_LIBS()
      FCFLAGS="$_gx_check_netcdf_save_FCFLAGS"
      unset _gx_check_netcdf_save_FCFLAGS
    ]
  )
])# _GX_CHECK_LIBS_NETCDF

# _GX_FORTRAN_PP_SRCEXT_PUSH(EXT) will ensure the Fortran test extension (stored
# in ac_ext) will cause the Fortran compiler to preprocess the test source file.
# Most Fortran compilers will preprocess the file based on the file extension,
# and of the known extension, F appears to work for all compilers.  This
# function accepts the file extension (without the preceeding .), similar
# to the AC_FC_SRCEXT and AC_FC_PP_SRCEXT macros.  If the extension is not
# provided, the default is 'F'.
# Unfortunately, this will not work for compilers that require a specific flag.
AC_DEFUN([_GX_FORTRAN_PP_SRCEXT_PUSH],[
  _gx_fortran_pp_srcext_save=$ac_ext
  AS_VAR_SET_IF([1], [ac_ext=$1], [ac_ext="F"])
])# _GX_FORTRAN_PP_SRCEXT_PUSH

# _GX_FORTRAN_PP_SRCEXT_POP() reverses the extension change done in
# _GX_FORTRAN_PP_SRCEXT_PUSH.  If this pop is called without a preceeding
# push, the default Fortran file extension (f) will be used.
AC_DEFUN([_GX_FORTRAN_PP_SRCEXT_POP], [
  AS_VAR_SET_IF([_gx_fortran_pp_srcext_save], [ac_ext=${_gx_fortran_pp_srcext_save}], [ac_ext="f"])
  unset _gx_fortran_pp_srcext_save
])# _GX_FORTRAN_PP_SRCEXT_POP

AC_DEFUN([_GX_CHECK_NETCDF_INC],[
  dnl First check if netcdf.inc is includable via the Fortran preprocessor
  dnl
  dnl Ensure the Fortran compiler will run the preprocessor
  _GX_FORTRAN_PP_SRCEXT_PUSH([F])
  AC_CACHE_CHECK([for netcdf.inc usability], [gx_cv_header_netcdf_inc], [
    AC_COMPILE_IFELSE(AC_LANG_PROGRAM([],
      [@%:@include <netcdf.inc>]),
      [gx_cv_header_netcdf_inc=yes],
      [AC_COMPILE_IFELSE(AC_LANG_PROGRAM([],
        [      include "netcdf.inc"],
        [gx_cv_header_netcdf_inc=yes],
        [gx_cv_header_netcdf_inc=no]))
    ])
  ])
  dnl Restore the SRCEXT
  _GX_FORTRAN_PP_SRCEXT_POP()
  AS_IF([test x"$gx_cv_header_netcdf_inc" = x"yes"],
    [AC_DEFINE([HAVE_NETCDF_INC], 1, [Define to 1 if netcdf.inc is available])])
])# _GX_CHECK_NETCDF_INC

AC_DEFUN([_GX_SEARCH_NETCDFF_LIBS], [
dnl First check for Fortran 77 bindings, even if AC_LANG is Fortran
AC_CACHE_CHECK([for nf_create in -lnetcdff], [gx_cv_search_nf_create],
[gx_search_netcdf_libs_save_LIBS=$LIBS
for gx_lib in '' 'netcdff'; do
  if test -z "$gx_lib"; then
    gx_res="none required"
  else
    gx_res=-l$gx_lib
    LIBS="$gx_res $gx_search_netcdf_libs_save_LIBS"
  fi
  AC_LINK_IFELSE(AC_LANG_PROGRAM([], [[      iret = nf_create("foo.nc", 1, ncid)]]),
    [gx_cv_search_nf_create="$gx_res"])
  AS_VAR_SET_IF([gx_cv_search_nf_create], [break])
done
AS_VAR_SET_IF([gx_cv_search_nf_create], , [gx_cv_search_nf_create=no])
AS_IF([test "$gx_cv_search_nf_create" != no],
  [test "$gx_cv_search_nf_create" = "none required" || LIBS="$ac_res $LIBS"])
])
dnl Check if Fortran 90 bindings are available in the netcdff library.  This
dnl check if similar to the one above.  If the library found above works,
dnl then this will add nothing to LIBS
AC_LANG_CASE(
[Fortran],[
AC_CACHE_CHECK([for nf90_create in -lnetcdff], [gx_cv_search_nf90_create],
[gx_search_netcdf_libs_save_LIBS=$LIBS
for gx_lib in '' 'netcdff'; do
  if test -z "$gx_lib"; then
    gx_res="none required"
  else
    gx_res=-l$gx_lib
    LIBS="$gx_res $gx_search_netcdf_libs_save_LIBS"
  fi
  AC_LINK_IFELSE(AC_LANG_PROGRAM([], [[      use netcdf
      iret = nf90_open("foo.nc", 1, ncid)]]),
    [gx_cv_search_nf90_create="$gx_res"])
  AS_VAR_SET_IF([gx_cv_search_nf90_create], [break])
done
AS_VAR_SET_IF([gx_cv_search_nf90_create], , [gx_cv_search_nf90_create=no])
AS_IF([test "$gx_cv_search_nf_create" != no],
  [test "gx_cv_search_nf90_create" = "none required" || LIBS="$ac_res $LIBS"])
])
])
])# _GX_SEARCH_NETCDFF_LIBS
