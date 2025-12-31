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

# Set common test settings.
. ../test-lib.sh

if [ -z "${parser_skip}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.yaml
title: test_prepend_date
base_date: 1 1 1 0 0 0
diag_files:
- file_name: test_non_static
  time_units: hours
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: average
    kind: r4
- file_name: test_static
  time_units: hours
  unlimdim: time
  freq: -1 hours
  varlist:
  - module: ocn_mod
    var_name: var2
    reduction: none
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc
my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n/" | cat > input.nml
test_expect_success "Running diag_manager and checking that the date was prepended correctly (test $my_test_count)" '
  mpirun -n 1 ../test_prepend_date
'

cat <<_EOF > diag_table.yaml
title: test_prepend_date
base_date: 1 1 1 0 0 0
diag_files:
- file_name: test_non_static
  time_units: hours
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
_EOF

printf "&diag_manager_nml \n use_modern_diag=.true. \n/ \n &test_prepend_date_nml \n pass_diag_time=.false. \n /" | cat > input.nml

test_expect_failure "Running diag_manager with fields that have a different start time (test $my_test_count)" '
  mpirun -n 1 ../test_prepend_date
'

fi
test_done
