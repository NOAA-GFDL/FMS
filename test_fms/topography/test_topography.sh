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
# execute tests in the test_fms/topography directory.

# Caitlyn McAllister

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
cat <<_EOF > input.nml
&topography_nml
  topog_file = 'topography.data.nc'
  water_file = 'water.data.nc'
  /
&gaussian_topog_nml
  height = 5000.0, 3000.0, 3000.0, 3000.0,
  olon   =   90.0,  255.0,  285.0,    0.0,
  olat   =   45.0,   45.0,  -15.0,  -90.0,
  wlon   =   15.0,   10.0,    5.0,  180.0,
  wlat   =   15.0,   25.0,   25.0,   20.0, 
  /
_EOF

touch input.nml

# Run the test.

test_expect_success "Test topography_mod: r4_kind" '
  mpirun -n 2 ./test_topography_r4
'

test_expect_success "Test topography_mod: r8_kind" '
  mpirun -n 2 ./test_topography_r8
'
test_done
