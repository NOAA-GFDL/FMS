# This is the main automake file for FMS.
# Ed Hartnett 2/21/2019

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = include platform constants tridiagonal mpp memutils fms \
affinity mosaic time_manager diag_manager drifters axis_utils horiz_interp \
time_interp column_diagnostics block_control \
data_override astronomy field_manager coupler diag_integral monin_obukhov \
interpolator fft amip_interp oda_tools exchange topography \
tracer_manager station_data sat_vapor_pres random_numbers libFMS \
test_fms/diag_manager test_fms/data_override test_fms/exchange \
test_fms/monin_obukhov test_fms/drifters test_fms/mosaic \
test_fms/interpolator test_fms/fms \
test_fms/mpp test_fms/mpp_io  test_fms/time_interp test_fms/time_manager \
test_fms/horiz_interp \
test_fms/field_manager test_fms/fft test_fms/axis_utils

EXTRA_DIST = README.md

CLEANFILES = *__genmod.mod *__genmod.f90
