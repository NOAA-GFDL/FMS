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

cat <<_EOF > diag_table.ens_01.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_ens
  time_units: days
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r8
_EOF

cat <<_EOF > diag_table.ens_02.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_ens
  time_units: days
  unlimdim: time
  freq: 2 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r8
_EOF

cat <<_EOF > input.nml
&diag_manager_nml
  use_modern_diag = .True.
/

&ensemble_nml
   ensemble_size = 2
/
_EOF

my_test_count=1
test_expect_success "Running diag_manager with 2 ensembles (test $my_test_count)" '
  mpirun -n 2 ../test_ens_runs
'

cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_ens
  time_units: days
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r8
_EOF

my_test_count=`expr $my_test_count + 1`
test_expect_failure "Running diag_manager with both diag_table.yaml and diag_table.ens_xx.yaml files present (test $my_test_count)" '
  mpirun -n 2 ../test_ens_runs
'

# Both ensembles have the same diag table yaml
cat <<_EOF > input.nml
&diag_manager_nml
  use_modern_diag = .True.
/

&ensemble_nml
   ensemble_size = 2
/
&test_ens_runs_nml
  using_ens_yaml = .false.
/
_EOF
rm -rf diag_table.ens_01.yaml diag_table.ens_02.yaml
my_test_count=`expr $my_test_count + 1`
test_expect_success "Running diag_manager with 2 ensembles, both ensembles have the same yaml (test $my_test_count)" '
  mpirun -n 2 ../test_ens_runs
'

fi
test_done
