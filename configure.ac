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


# These files will be created when the configure script is run.
AC_CONFIG_FILES([Makefile
                 ])
AC_OUTPUT()
