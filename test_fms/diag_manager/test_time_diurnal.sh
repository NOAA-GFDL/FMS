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
# tests the diurnal (daily average) reduction method

# Set common test settings.
. ../test-lib.sh

if [ -z "${skipflag}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.yaml
title: test_diurnal
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_diurnal
  time_units: hours
  unlimdim: time
  freq: 1 months
  varlist:
  - module: ocn_mod
    var_name: var4
    output_name: var4
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var2
    output_name: var2
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var1
    output_name: var1
    reduction: diurnal3
    kind: r4
- file_name: test_diurnal_regional
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: latlon
    corner1: 78. 78.
    corner2: 78. 78.
    corner3: 81. 81.
    corner4: 81. 81.
  freq: 1 months
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_diurnal
    reduction: diurnal3
    zbounds: 2. 3.
    kind: r4
_EOF

export OMP_NUM_THREADS=1

printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n" > input.nml

test_expect_success "monthly simple diurnal output" '
  mpirun -n 6 ../test_diag_diurnal
'

test_expect_success "checking results for diurnal test simple" '
  mpirun -n 1 ../check_time_diurnal
'

printf "&test_diag_diurnal_nml \n test_case=0 \n mask_case=1 \n / \n" >> input.nml

test_expect_success "monthly diurnal output with logical mask" '
  mpirun -n 6 ../test_diag_diurnal
'
test_expect_success "checking results for diurnal test with logical mask" '
  mpirun -n 1 ../check_time_diurnal
'

printf "&test_diag_diurnal_nml \n test_case=0 \n mask_case=2 \n / \n" >> input.nml

test_expect_success "monthly diurnal output with real mask" '
  mpirun -n 6 ../test_diag_diurnal
'
test_expect_success "checking results for diurnal test with real mask" '
  mpirun -n 1 ../check_time_diurnal
'

export OMP_NUM_THREADS=2

printf "&test_diag_diurnal_nml \n test_case=1 \n / \n" >> input.nml

test_expect_success "monthly diurnal output with openmp" '
  mpirun -n 6 ../test_diag_diurnal
'
test_expect_success "checking results for diurnal test with openmp" '
  mpirun -n 1 ../check_time_diurnal
'

printf "&test_diag_diurnal_nml \n test_case=1 \n mask_case=1 \n / \n" >> input.nml

test_expect_success "monthly diurnal output with openmp and real mask" '
  mpirun -n 6 ../test_diag_diurnal
'
test_expect_success "checking results for diurnal test with openmp and real mask" '
  mpirun -n 1 ../check_time_diurnal
'

printf "&test_diag_diurnal_nml \n test_case=1 \n mask_case=2 \n / \n" >> input.nml

test_expect_success "monthly diurnal output with openmp and logical mask" '
  mpirun -n 6 ../test_diag_diurnal
'
test_expect_success "checking results for diurnal test with openmp and logical mask" '
  mpirun -n 1 ../check_time_diurnal
'

fi

test_done
