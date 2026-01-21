#!/bin/sh

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

# Copyright (c) 2019-2020 Ed Hartnett, Seth Underwood

# Set common test settings.
. ../test-lib.sh

if [ -z "${parser_skip}" ]; then
# create and enter directory for in/output files
output_dir

#TODO replace with yaml diag_table and set diag_manager_nml::use_modern_diag=.true.
cat <<_EOF > diag_table.yaml
title: test_max
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_max
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    output_name: var0_max
    reduction: max
    kind: r4
  - module: ocn_mod
    var_name: var1
    output_name: var1_max
    reduction: max
    kind: r4
  - module: ocn_mod
    var_name: var2
    output_name: var2_max
    reduction: max
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3_max
    reduction: max
    kind: r4
  - module: ocn_z_mod
    var_name: var4
    output_name: var4_max
    reduction: max
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3_Z_max
    reduction: max
    zbounds: 2. 3.
    kind: r4
- file_name: test_max_regional
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: latlon
    corner1: 78. 78.
    corner2: 78. 78.
    corner3: 81. 81.
    corner4: 81. 81.
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_max
    reduction: max
    zbounds: 2. 3.
    kind: r4
_EOF

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 0 \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n &test_reduction_methods_nml \n test_case = 0 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n &test_reduction_methods_nml \n test_case = 0 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'

export OMP_NUM_THREADS=2
my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method with openmp (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method with openmp (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'
export OMP_NUM_THREADS=1

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method with halo output (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method with halo output (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "max" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "max" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_max
'
fi
test_done