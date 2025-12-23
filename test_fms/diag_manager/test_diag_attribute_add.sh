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
- file_name: food_file
  freq: 4 hours
  time_units: hours
  unlimdim: time
  reduction: none
  kind: r4
  module: food_mod
  varlist:
    - var_name: potatoes
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc

my_test_count=1
cat <<_EOF > input.nml
&diag_manager_nml
 use_modern_diag=.true.
 max_field_attributes = 10
/
_EOF

test_expect_success "Testing diag_field_attribute_add (test $my_test_count)" '
  mpirun -n 1 ../test_diag_attribute_add
'
fi
test_done
