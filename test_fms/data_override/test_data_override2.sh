#!/bin/sh

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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/data_override directory.

# Ed Hartnett 11/26/19
# Uriel Ramirez 07/22/20

# Set common test settings.
. ../test_common.sh

# Run the ongrid test case with 2 halos in x and y
touch input.nml
printf '"OCN", "runoff", "runoff", "./INPUT/runoff.daitren.clim.1440x1080.v20180328.nc", "none" ,  1.0' | cat > data_table
[ ! -d "INPUT" ] && mkdir -p "INPUT"
run_test test_data_override_ongrid 6
rm -rf "INPUT"

# Run the ongrid test case again with no halos
printf "&test_data_override_ongrid_nml \n nhalox=0 \n nhaloy=0\n/" | cat > input.nml
[ ! -d "INPUT" ] && mkdir -p "INPUT"
run_test test_data_override_ongrid 6
rm -rf "INPUT"

# Copy to builddir and rename data files for tests.
#cp $top_srcdir/test_fms/data_override/data_table_base data_table
#cp $top_srcdir/test_fms/data_override/diag_table_base diag_table
#cp -r $top_srcdir/test_fms/data_override/INPUT $top_builddir/test_fms/data_override/INPUT

# Both tests are skipped in the bats file, so commented out here.
#tnum=$( printf "%2.2d" 1 )
#sed "s/<test_num>/${tnum}/"  $top_srcdir/test_fms/data_override/input_base.nml > input.nml
#run_test test_data_override 2 skip

#tnum=$( printf "%2.2d" 2 )
#sed "s/<test_num>/${tnum}/"  $top_srcdir/test_fms/data_override/input_base.nml > input.nml
#run_test test_data_override 2 skip

