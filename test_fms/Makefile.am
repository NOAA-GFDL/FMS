#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
#***********************************************************************

# This is the automake file for the test_fms directory.
# Ed Hartnett 9/20/2019

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# Make targets will be run in each subdirectory. Order is significant.
SUBDIRS = astronomy coupler diag_manager data_override exchange monin_obukhov drifters \
mosaic2 interpolator fms mpp time_interp time_manager horiz_interp topography \
field_manager axis_utils affinity fms2_io parser string_utils sat_vapor_pres tracer_manager \
random_numbers diag_integral column_diagnostics tridiagonal block_control

# testing utility scripts to distribute
EXTRA_DIST = test-lib.sh.in intel_coverage.sh.in tap-driver.sh
