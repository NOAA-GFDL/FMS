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
# execute tests in the test_fms/coupler directory.

# Ed Hartnett 11/26/19
# Uriel Ramirez 07/22/20

# Set common test settings.
. ../test-lib.sh

rm -f input.nml
touch input.nml

# diag_table for test
cat <<_EOF > diag_table
test_coupler
1 1 1 0 0 0
#output files
 "coupler_types_bc2",  1, "days", 1, "days", "time"
 "coupler_types_bc1",  1, "days", 1, "days", "time"
#output variables
 "test_coupler_types", "bc1_var2d_1", "bc1_variable_2d_1_min", "coupler_types_bc1", "all", "min", "none", 2
 "test_coupler_types", "bc1_var2d_2", "bc1_variable_2d_2_max", "coupler_types_bc1", "all", "max", "none", 2
 "test_coupler_types", "bc1_var3d_1", "bc1_variable_3d_1", "coupler_types_bc1", "all", "rms", "none", 2
 "test_coupler_types", "bc1_var3d_2", "bc1_variable_3d_2", "coupler_types_bc1", "all", "avg", "none", 2
 "test_coupler_types", "bc2_var2d_1", "bc2_variable_2d_1_min", "coupler_types_bc2", "all", "min", "none", 2
 "test_coupler_types", "bc2_var2d_2", "bc2_variable_2d_2_max", "coupler_types_bc2", "all", "max", "none", 2
 "test_coupler_types", "bc2_var3d_1", "bc2_variable_3d_1", "coupler_types_bc2", "all", "rms", "none", 2
 "test_coupler_types", "bc2_var3d_2", "bc2_variable_3d_2", "coupler_types_bc2", "all", "avg", "none", 2
_EOF
# we'll just make both in case compiled with yaml support
cat <<_EOF > diag_table.yaml
title: test_coupler
base_date: 1 1 1 0 0 0
diag_files:
- file_name: coupler_types_bc2
  filename_time: end
  freq: 1 days
  time_units: days
  unlimdim: time
  varlist:
  - module: test_coupler_types
    var_name: bc1_var2d_1
    output_name: bc1_variable_2d_1_min
    reduction: min
  - module: test_coupler_types
    var_name: bc1_var2d_2
    output_name: bc1_variable_2d_2_max
    reduction: max
  - module: test_coupler_types
    var_name: bc1_var3d_1
    output_name: bc1_variable_3d_1
    reduction: rms
  - module: test_coupler_types
    var_name: bc1_var3d_2
    output_name: bc1_variable_3d_2
    reduction: avg
- file_name: coupler_types_bc1
  filename_time: end
  freq: 1 days
  time_units: days
  unlimdim: time
  varlist:
  - module: test_coupler_types
    var_name: bc2_var2d_1
    output_name: bc2_variable_2d_1_min
    reduction: min
  - module: test_coupler_types
    var_name: bc2_var2d_2
    output_name: bc2_variable_2d_2_max
    reduction: max
  - module: test_coupler_types
    var_name: bc2_var3d_1
    output_name: bc2_variable_3d_1
    reduction: rms
  - module: test_coupler_types
    var_name: bc2_var3d_2
    output_name: bc2_variable_3d_2
    reduction: avg
_EOF

cat <<_EOF > data_table
"ATM", "bc1_var2d_1",  "bc1_variable_2d_1_min", "coupler_types_bc1.nc", .false., 300.0
_EOF

rm -rf INPUT
mkdir INPUT


test_expect_success "coupler types interfaces (r4_kind)" '
  mpirun -n 4 ./test_coupler_types_r4
'

test_expect_success "coupler types interfaces (r8_kind)" '
  mpirun -n 4 ./test_coupler_types_r8
'

# delete lines from the table to make sure we see the difference in the send_data return status
sed -i '8,12{d}' diag_table
sed -i '10,13{d}' diag_table.yaml
sed -i '18,25{d}' diag_table.yaml
cat <<_EOF > input.nml
&test_coupler_types_nml
   fail_return_status=.true.
/
_EOF


test_expect_success "coupler types interfaces - check send_data return vals (r4_kind)" '
  mpirun -n 4 ./test_coupler_types_r4
'

test_expect_success "coupler types interfaces - check send_data return vals (r8_kind)" '
  mpirun -n 4 ./test_coupler_types_r8
'

mkdir RESTART

test_expect_success "coupler register restart 2D(r4_kind)" '
  mpirun -n 1 ./test_coupler_2d_r4
'
test_expect_success "coupler register restart 2D(r8_kind)" '
  mpirun -n 1 ./test_coupler_2d_r8
'

test_expect_success "coupler register restart 3D (r4_kind)" '
  mpirun -n 1 ./test_coupler_3d_r4
'

test_expect_success "coupler register restart 3D (r8_kind)" '
  mpirun -n 1 ./test_coupler_3d_r8
'

test_expect_success "test atmos_ocean_fluxes (r4_kind)" '
  mpirun -n 1 ./test_atmos_ocean_fluxes_r4
'

test_expect_success "test atmos_ocean_fluxes (r8_kind)" '
  mpirun -n 1 ./test_atmos_ocean_fluxes_r8
'

rm -rf RESTART
test_done
