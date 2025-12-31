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
[ ! -d "INPUT" ] && mkdir -p "INPUT"

cat <<_EOF > data_table.ens_01.yaml
data_table:
 - grid_name: OCN
   fieldname_in_model: runoff
   override_file:
   - fieldname_in_file: runoff
     file_name: INPUT/runoff.daitren.clim.1440x1080.v20180328_ens_01.nc
     interp_method: none
   factor: 1.0
_EOF

cat <<_EOF > data_table.ens_02.yaml
data_table:
 - grid_name: OCN
   fieldname_in_model: runoff
   override_file:
   - fieldname_in_file: runoff
     file_name: INPUT/runoff.daitren.clim.1440x1080.v20180328_ens_02.nc
     interp_method: none
   factor: 1.0
_EOF

cat <<_EOF > input_base.nml
&data_override_nml
  use_data_table_yaml = .True.
/

&test_data_override_ongrid_nml
  test_case = 5
  write_only = .False.
/

&ensemble_nml
   ensemble_size = 2
/
_EOF

#The test only runs with yaml
if [ -z $parser_skip ]; then
  rm -rf INPUT/.

  sed 's/write_only = .False./write_only = .True./g' input_base.nml > input.nml
  test_expect_success "Creating input files" '
    mpirun -n 12 ../test_data_override_ongrid
  '

  cp input_base.nml input.nml
  test_expect_success "test_data_override with two ensembles  -yaml" '
    mpirun -n 12 ../test_data_override_ongrid
  '

  cat <<_EOF > data_table.yaml
data_table:
 - grid_name: OCN
   fieldname_in_model: runoff
   override_file:
   - fieldname_in_file: runoff
     file_name: INPUT/runoff.daitren.clim.1440x1080.v20180328_ens_01.nc
     interp_method: none
   factor: 1.0
_EOF

  test_expect_failure "test_data_override with both data_table.yaml and data_table.ens_xx.yaml files" '
    mpirun -n 12 ../test_data_override_ongrid
  '

  cat <<_EOF > input.nml
&data_override_nml
  use_data_table_yaml = .True.
/

&test_data_override_ongrid_nml
  test_case = 6
  write_only = .False.
/

&ensemble_nml
  ensemble_size = 2
/
_EOF

  rm -rf INPUT/.
  rm -rf data_table.ens_01.yaml data_table.ens_02.yaml

  test_expect_success "test_data_override with two ensembles, same yaml file" '
    mpirun -n 12 ../test_data_override_ongrid
  '

  rm -rf INPUT
fi

test_done
