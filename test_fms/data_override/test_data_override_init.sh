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
#
# Copyright (c) 2019-2021 Ed Hartnett, Uriel Ramirez, Seth Underwood

# Set common test settings.
. ../test-lib.sh

output_dir

# data_override with the default table (not setting namelist)
rm -rf input.nml data_table data_table.yaml
touch input.nml
touch data_table
test_expect_success "data_override_init with the default table" '
 mpirun -n 1 ../test_data_override_init
'

cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.false.
/
_EOF
test_expect_success "data_override_init setting use_data_table_yaml = .false." '
 mpirun -n 1 ../test_data_override_init
'

touch data_table.yaml
test_expect_failure "data_override_init both tables present" '
 mpirun -n 1 ../test_data_override_init
'
if [ ! -z $parser_skip ]; then
rm -rf data_table.yaml
cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.true.
/
_EOF
  test_expect_failure "data_override_init setting use_data_table_yaml = .true. but no compiling with yaml" '
    mpirun -n 1 ../test_data_override_init
  '
else
rm -rf data_table
cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.true.
/
_EOF
test_expect_success "data_override_init setting use_data_table_yaml = .true." '
 mpirun -n 1 ../test_data_override_init
'
fi
test_done