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
rm -rf data_table data_table.yaml input.nml input_base.nml

if [ -z $parser_skip ]; then
  use_yaml=true
  cat <<_EOF > data_table.yaml
data_table:
 - grid_name: OCN
   fieldname_in_model: runoff
   override_file:
   - fieldname_in_file: runoff
     file_name: INPUT/runoff.daitren.clim.1440x1080.v20180328.nc
     interp_method: none
   factor: 1.0
_EOF
else
  use_yaml=false
  printf '"OCN", "runoff", "runoff", "./INPUT/runoff.daitren.clim.1440x1080.v20180328.nc", "none" ,  1.0' | cat > data_table
fi

prepare_input_nml () {
  halo_size=$1
  write_only=$2
  init_with_mode=$3

  if [ -z $init_with_mode ]
  then
    init_with_mode=false
  fi

  cat <<_EOF > input.nml
&data_override_nml
  use_data_table_yaml=.${use_yaml}.
/

&test_data_override_ongrid_nml
  nhalox=${halo_size}
  nhaloy=${halo_size}
  write_only = .${write_only}.
  init_with_mode=.${init_with_mode}.
/
_EOF
}

[ ! -d "INPUT" ] && mkdir -p "INPUT"

prepare_input_nml 2 true
test_expect_success "Creating input files" '
  mpirun -n 6 ../test_data_override_ongrid
'

prepare_input_nml 2 false false
test_expect_success "data_override on grid with 2 halos in x and y" '
  mpirun -n 6 ../test_data_override_ongrid
'

prepare_input_nml 2 false true
test_expect_success "data_override on grid with 2 halos in x and y (init with explicit modes)" '
  mpirun -n 6 ../test_data_override_ongrid
'

prepare_input_nml 0 false false
test_expect_success "data_override on grid with 0 halos in x and y" '
  mpirun -n 6 ../test_data_override_ongrid
'

prepare_input_nml 0 false true
test_expect_success "data_override on grid with 0 halos in x and y (init with explicit modes)" '
  mpirun -n 6 ../test_data_override_ongrid
'

for KIND in r4 r8
do
  rm -rf INPUT/*

  test_expect_success "data_override get_grid_v1 (${KIND})" '
    mpirun -n 1 ../test_get_grid_v1_${KIND}
  '
done

rm -rf INPUT *.nc # remove any leftover files to reduce size

test_done
