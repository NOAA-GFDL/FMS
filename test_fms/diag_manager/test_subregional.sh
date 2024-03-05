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

# Set common test settings.
. ../test-lib.sh

if [ -z "${skipflag}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.yaml
title: test_subregional
base_date: 2 1 1 0 0 0

diag_files:
# This is to test a file with multiple z axis
- file_name: test_subZaxis
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_Z1
    reduction: none
    kind: r4
    zbounds: 2. 3.
  - module: ocn_mod
    var_name: var3
    output_name: var3_Z2
    reduction: none
    kind: r4
    zbounds: 3. 5.
- file_name: test_subregional
  freq: 6 hours
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: latlon
    corner1: 60. 60.
    corner2: 60. 65.
    corner3: 65. 65.
    corner4: 65. 60.
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_min
    reduction: min
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3_max
    reduction: max
    kind: r4
- file_name: test_subregional2
  freq: 6 hours
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: index
    corner1: 60 60
    corner2: 60 65
    corner3: 65 65
    corner4: 65 60
    tile: 1
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_min
    reduction: min
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3_max
    reduction: max
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n/" | cat > input.nml
test_expect_success "Running diag_manager with different subregions (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'

cat <<_EOF > diag_table.yaml
title: test_corner_subregional
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_corner1
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var2c
    output_name: var2c_avg
    reduction: average
    kind: r4
  sub_region:
  - grid_type: latlon
    corner1: 17. 17.
    corner2: 17. 20.
    corner3: 20. 17.
    corner4: 20. 20.
- file_name: test_corner2
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var2c
    output_name: var2c_avg
    reduction: average
    kind: r4
  sub_region:
  - grid_type: latlon
    corner1: 17. 17.
    corner2: 20. 17.
    corner3: 17. 17.
    corner4: 20. 17.
- file_name: test_corner3
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var2c
    output_name: var2c_avg
    reduction: average
    kind: r4
  sub_region:
  - grid_type: latlon
    corner1: 17. 17.
    corner2: 20. 17.
    corner3: 17. 33.
    corner4: 20. 33.
_EOF

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n/" | cat > input.nml
test_expect_success "Running diag_manager with corner diagnotics (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'

my_test_count=`expr $my_test_count + 1`
test_expect_success "Checking results from diag_manager with different subregions (test $my_test_count)" '
  mpirun -n 1 ../check_subregional
'

fi
test_done
