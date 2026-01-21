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
SKIP_TESTS="test_data_override_multi_file.2"

output_dir
rm -rf data_table data_table.yaml input.nml input_base.nml

if [ -z $parser_skip ]; then

cat <<_EOF > input_base.nml
&data_override_nml
use_data_table_yaml=.True.
/
&test_data_override_ongrid_nml
  test_case = 7
  nlon = 5
  nlat = 6
  layout = 1, 2
  write_only = .False.
/
_EOF

cat <<_EOF > data_table.yaml
data_table:
 - grid_name           : OCN
   fieldname_in_model  : runoff
   override_file:
   - file_name         : INPUT/hadisst_ice.data_yr1.nc
     fieldname_in_file : runoff
     interp_method     : none
     multi_file:
     - next_file_name: INPUT/hadisst_ice.data_yr2.nc
       prev_file_name: INPUT/hadisst_ice.data_yr0.nc
   factor              : 1
_EOF

[ ! -d "INPUT" ] && mkdir -p "INPUT"
rm -rf INPUT/*

sed 's/write_only = .False./write_only = .True./g' input_base.nml > input.nml
test_expect_success "Creating input files" '
  mpirun -n 2 ../test_data_override_ongrid
'

cp input_base.nml input.nml
test_expect_success "data_override multi_file" '
  mpirun -n 2 ../test_data_override_ongrid
'

rm -rf INPUT *.nc # remove any leftover files to reduce size

fi

test_done
