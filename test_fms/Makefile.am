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

# This is the automake file for the test_fms directory.
# Ed Hartnett 9/20/2019

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = \
  affinity \
  axis_utils \
  data_override \
  diag_manager \
  drifters \
  exchange
#  field_manager \
#  fms \
#  fms2_io \
#  horiz_interp \
#  interpolator \
#  monin_obukhov \
#  mosaic \
#  mpp \
#  mpp_io \
#  time_interp \
#  time_manager



# This input file must be distributed, it is turned into
# test_common.sh by configure.
EXTRA_DIST = test-lib.sh.in
