# This is the main automake file for FMS.
# Ed Hartnett 2/21/2019

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = include platform constants tridiagonal mpp memutils fms \
affinity fms2_io mosaic time_manager axis_utils diag_manager drifters \
horiz_interp time_interp column_diagnostics block_control \
data_override astronomy field_manager coupler diag_integral \
monin_obukhov interpolator amip_interp exchange topography \
tracer_manager sat_vapor_pres random_numbers libFMS test_fms

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = FMS.pc

EXTRA_DIST = README.md

CLEANFILES = *__genmod.mod *__genmod.f90
