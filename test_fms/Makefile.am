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
SUBDIRS = diag_manager data_override exchange monin_obukhov drifters \
mosaic interpolator fms fms2_io mpp mpp_io time_interp time_manager \
horiz_interp field_manager axis_utils affinity

# This input file must be distributed, it is turned into
# test_common.sh by configure.
EXTRA_DIST = test_common.sh.in
