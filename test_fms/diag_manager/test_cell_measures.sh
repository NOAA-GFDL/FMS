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
title: test_diag_manager
base_date: 2 1 1 0 0 0

diag_files:
# Here file 1 does not have the "area" variable so the associated files attribute is expected
- file_name: file1
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: fun_mod
    var_name: var1
    reduction: average
    kind: r4
- file_name: file2
  freq: 1 hours
  time_units: hours
  unlimdim: time
  module: fun_mod
  reduction: none
  kind: r4
  varlist:
  - var_name: var1
  - var_name: area
    output_name: area_file2
- file_name: static_file
  freq: -1
  time_units: hours
  unlimdim: time
  varlist:
  - module: fun_mod
    var_name: area
    reduction: none
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n/" | cat > input.nml
test_expect_success "Running diag_manager with fields with cell measures (area, volume) (test $my_test_count)" '
  mpirun -n 1 ../test_cell_measures
'
fi
test_done
