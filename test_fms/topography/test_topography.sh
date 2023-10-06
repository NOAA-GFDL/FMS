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
# execute tests in the test_fms/astronomy directory.

# Caitlyn McAllister

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
cat <<EOF > input.nml
&topography_nml
  topog_file = "topography.data.nc",
  water_file = "water.data.nc"
/

EOF

# Run the test.

test_expect_success "Test topography: r4_kind" '
  mpirun -n 1 ./test_topography_r4
'

sync; rm -f *.nc

test_expect_success "Test topography: r8_kind" '
  mpirun -n 1 ./test_topography_r8
'

rm -f *.nc

test_done
