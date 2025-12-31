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

cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_0days
  time_units: days
  unlimdim: time
  freq: 0 days
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r4
  - module: ocn_mod
    var_name: var2
    reduction: none
    kind: r4
  - module: ocn_mod
    var_name: var0_openmp
    reduction: none
    kind: r4
_EOF

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 0 \n/" | cat > input.nml
test_expect_success "Running diag_manager with 0 days frequency (test $my_test_count)" '
  mpirun -n 1 ../test_output_every_freq
'

cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_0days
  time_units: days
  unlimdim: time
  freq: 0 days
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r4
  - module: ocn_mod
    var_name: var1
    reduction: none
    kind: r4
  - module: ocn_mod
    var_name: var2
    reduction: none
    kind: r4
_EOF

my_test_count=`expr $my_test_count + 1`
test_expect_failure "Running diag_manager with 0 days frequency and variable are calling send data at different frequencies (test $my_test_count)" '
  mpirun -n 1 ../test_output_every_freq
'

cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_0days
  time_units: days
  unlimdim: time
  freq: 0 days
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: average
    kind: r4
_EOF

my_test_count=`expr $my_test_count + 1`
test_expect_failure "Running diag_manager with 0 days frequency but using average reduction method (test $my_test_count)" '
  mpirun -n 1 ../test_output_every_freq
'

cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_0days
  time_units: days
  unlimdim: time
  freq: -1 days
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: average
    kind: r4
_EOF

my_test_count=`expr $my_test_count + 1`
test_expect_failure "Running diag_manager with -1 days frequency but using average reduction method (test $my_test_count)" '
  mpirun -n 1 ../test_output_every_freq
'
fi
test_done
