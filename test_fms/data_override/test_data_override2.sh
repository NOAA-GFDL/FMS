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
#
# Copyright (c) 2019-2021 Ed Hartnett, Uriel Ramirez, Seth Underwood

# Set common test settings.
. ../test-lib.sh

# Tests to skip
SKIP_TESTS="$(basename $0 .sh).[1-2]"

setup_test_dir () {
  local halo_size
  test "$#" = 1 && { halo_size=$1; } ||
  BUG "required parameter for halo size not present"
  rm -rf data_table input.nml INPUT
  cat <<_EOF > data_table
"OCN", "runoff", "runoff", "./INPUT/runoff.daitren.clim.1440x1080.v20180328.nc", "none" ,  1.0
_EOF

cat <<_EOF > input.nml
&test_data_override_ongrid_nml
  nhalox=${halo_size}
  nhaloy=${halo_size}
/
_EOF
}

setup_test_dir 2
test_expect_success "data_override on grid with 2 halos in x and y" '
  mpirun -n 6 ./test_data_override_ongrid
'

setup_test_dir 0
test_expect_success "data_override on grid with no halos" '
  mpirun -n 6 ./test_data_override_ongrid
'

# Run the get_grid_v1 test:
test_expect_success "data_override get_grid_v1" '
  mpirun -n 1 ./test_get_grid_v1
'

test_done
