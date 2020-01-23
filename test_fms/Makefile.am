# This is the automake file for the test_fms directory.
# Ed Hartnett 9/20/2019

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = diag_manager data_override exchange monin_obukhov drifters \
mosaic interpolator fms fms2_io mpp mpp_io time_interp time_manager \
horiz_interp field_manager axis_utils affinity

# This input file must be distributed, it is turned into
# test_common.sh by configure.
EXTRA_DIST = test_common.sh.in