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

# TODO: Enable this test once generalized indices work is complete
SKIP_TESTS="test_data_override2_mono.2"

output_dir
[ ! -d "INPUT" ] && mkdir -p "INPUT"

cat <<_EOF > input_base.nml
&test_data_override_ongrid_nml
  test_case = 2
  write_only = .False.
/
_EOF

cat <<_EOF > data_table
"OCN", "runoff_increasing", "runoff", "./INPUT/bilinear_increasing.nc", "bilinear" ,  1.0
"OCN", "runoff_decreasing", "runoff", "./INPUT/bilinear_decreasing.nc", "bilinear" ,  1.0
_EOF

rm -rf INPUT/*
sed 's/write_only = .False./write_only = .True./g' input_base.nml > input.nml
test_expect_success "Creating input files" '
  mpirun -n 6 ../test_data_override_ongrid
'

cp input_base.nml input.nml
test_expect_success "test_data_override with monotonically increasing and decreasing data sets" '
  mpirun -n 6 ../test_data_override_ongrid
'
rm -f INPUT/*
sync

rm -rf data_table

cat <<_EOF > input_base.nml
&test_data_override_ongrid_nml
  test_case = 2
  write_only = .False.
/
&data_override_nml
  use_data_table_yaml = .True.
/
_EOF

cat <<_EOF > data_table.yaml
data_table:
- grid_name: OCN
  fieldname_in_model: runoff_increasing
  override_file:
  - fieldname_in_file: runoff
    file_name: ./INPUT/bilinear_increasing.nc
    interp_method: bilinear
  factor: 1.0
- grid_name: OCN
  fieldname_in_model: runoff_decreasing
  override_file:
  - fieldname_in_file: runoff
    file_name: ./INPUT/bilinear_decreasing.nc
    interp_method: bilinear
  factor: 1.0
_EOF

#Repeat the test with yaml if needed
if [ -z $parser_skip ]; then
  # TODO: Enable this test once generalized indices work is complete
  SKIP_TESTS="$SKIP_TESTS test_data_override2_mono.4"

  rm -rf INPUT/*
  sed 's/write_only = .False./write_only = .True./g' input_base.nml > input.nml
  test_expect_success "Creating input files" '
    mpirun -n 6 ../test_data_override_ongrid
  '

  cp input_base.nml input.nml
  test_expect_success "test_data_override with monotonically increasing and decreasing data sets  -yaml" '
    mpirun -n 6 ../test_data_override_ongrid
  '
fi

rm -rf INPUT
test_done
