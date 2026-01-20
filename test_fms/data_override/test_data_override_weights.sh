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
SKIP_TESTS="test_data_override_weights.2"

output_dir
[ ! -d "INPUT" ] && mkdir -p "INPUT"

cat <<_EOF > data_table.yaml
data_table:
- grid_name: OCN
  fieldname_in_model: runoff_obs
  override_file:
  - fieldname_in_file: runoff
    file_name: ./INPUT/bilinear_increasing.nc
    interp_method: bilinear
  factor: 1.0
- grid_name: OCN
  fieldname_in_model: runoff_obs_weights
  override_file:
  - fieldname_in_file: runoff
    file_name: ./INPUT/bilinear_increasing.nc
    interp_method: bilinear
    external_weights:
    - file_name: ./INPUT/remap_file.nc
      source: fregrid
  factor: 1.0
_EOF

cat <<_EOF > input_base.nml
&data_override_nml
  use_data_table_yaml = .True.
/

&test_data_override_ongrid_nml
  test_case = 4
  nlon = 5
  nlat = 6
  layout = 1, 2
  write_only = .False.
/
_EOF

#The test only runs with yaml
if [ -z $parser_skip ]; then
  rm -rf INPUT/.

  sed 's/write_only = .False./write_only = .True./g' input_base.nml > input.nml
  test_expect_success "Creating input files" '
    mpirun -n 2 ../test_data_override_ongrid
  '

  cp input_base.nml input.nml
  test_expect_success "test_data_override with and without weight files  -yaml" '
    mpirun -n 2 ../test_data_override_ongrid
  '
fi

rm -rf INPUT
test_done
