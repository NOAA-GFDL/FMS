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

if [ -z "${parser_skip}" ] && [ -z "${parallel_skip}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.yaml
title: test_none
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_none
  freq: 6 hours
  time_units: hours
  unlimdim: time
  module: ocn_mod
  reduction: none
  kind: r4
  use_collective_writes: True
  varlist:
  - var_name: var0
    output_name: var0_none
  - var_name: var1
    output_name: var1_none
  - var_name: var2
    output_name: var2_none
  - var_name: var3
    output_name: var3_none
  - var_name: var4
    output_name: var4_none
    module: ocn_z_mod
  - var_name: var3
    output_name: var3_Z
    zbounds: 2. 3.
- file_name: test_none_regional
  freq: 6 hours
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: latlon
    corner1: 78. 78.
    corner2: 78. 78.
    corner3: 81. 81.
    corner4: 81. 81.
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_none
    reduction: none
    zbounds: 2. 3.
    kind: r4
_EOF

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 0 \n/" | cat > input.nml
test_expect_success "Running diag_manager with "none" reduction method - collective netcdf writes (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "none" reduction method (test $my_test_count)" '
  mpirun -n 1 ../check_time_none
'

if [ -z "${SKIP_NCDUMP_CHECKS}" ]; then
test_expect_success "Checking chunksizes" '
    ncdump -hs test_none.nc | grep -E "var4_none:_ChunkSizes *= *1, *2, *5, *16, *96" &&
    ncdump -hs test_none.nc | grep -E "var3_none:_ChunkSizes *= *1, *5, *16, *96" &&
    ncdump -hs test_none.nc | grep -E "var2_none:_ChunkSizes *= *1, *16, *96"
  '
fi

# Repeat the test with specified chunksizes
cat <<_EOF > diag_table.yaml
title: test_none
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_none
  freq: 6 hours
  time_units: hours
  unlimdim: time
  module: ocn_mod
  reduction: none
  kind: r4
  use_collective_writes: True
  chunksizes: 96, 8, 1, 1, 1
  varlist:
  - var_name: var0
    output_name: var0_none
  - var_name: var1
    output_name: var1_none
  - var_name: var2
    output_name: var2_none
  - var_name: var3
    output_name: var3_none
  - var_name: var4
    output_name: var4_none
    module: ocn_z_mod
  - var_name: var3
    output_name: var3_Z
    zbounds: 2. 3.
- file_name: test_none_regional
  freq: 6 hours
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: latlon
    corner1: 78. 78.
    corner2: 78. 78.
    corner3: 81. 81.
    corner4: 81. 81.
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_none
    reduction: none
    zbounds: 2. 3.
    kind: r4
_EOF

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 0 \n/" | cat > input.nml
test_expect_success "Running diag_manager with "none" reduction method - collective netcdf writes - specified chunksizes(test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "none" reduction method (test $my_test_count)" '
  mpirun -n 1 ../check_time_none
'

if [ -z "${SKIP_NCDUMP_CHECKS}" ]; then
test_expect_success "Checking chunksizes" '
    ncdump -hs test_none.nc | grep -E "var4_none:_ChunkSizes *= *1, *1, *1, *8, *96" &&
    ncdump -hs test_none.nc | grep -E "var3_none:_ChunkSizes *= *1, *1, *8, *96" &&
    ncdump -hs test_none.nc | grep -E "var2_none:_ChunkSizes *= *1, *8, *96"
  '
fi

# Repeat the test with specified chunksizes at variable level
cat <<_EOF > diag_table.yaml
title: test_none
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_none
  freq: 6 hours
  time_units: hours
  unlimdim: time
  module: ocn_mod
  reduction: none
  kind: r4
  use_collective_writes: True
  chunksizes: 96, 8, 1, 1, 1
  varlist:
  - var_name: var0
    output_name: var0_none
  - var_name: var1
    output_name: var1_none
  - var_name: var2
    output_name: var2_none
  - var_name: var3
    output_name: var3_none
  - var_name: var4
    output_name: var4_none
    module: ocn_z_mod
    chunksizes: 96, 8, 5, 2, 1
  - var_name: var3
    output_name: var3_Z
    zbounds: 2. 3.
- file_name: test_none_regional
  freq: 6 hours
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: latlon
    corner1: 78. 78.
    corner2: 78. 78.
    corner3: 81. 81.
    corner4: 81. 81.
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_none
    reduction: none
    zbounds: 2. 3.
    kind: r4
_EOF

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 0 \n/" | cat > input.nml
test_expect_success "Running diag_manager with "none" reduction method - collective netcdf writes - specified chunksizes at variable level (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "none" reduction method (test $my_test_count)" '
  mpirun -n 1 ../check_time_none
'

if [ -z "${SKIP_NCDUMP_CHECKS}" ]; then
test_expect_success "Checking chunksizes" '
    ncdump -hs test_none.nc | grep -E "var4_none:_ChunkSizes *= *1, *2, *5, *8, *96" &&
    ncdump -hs test_none.nc | grep -E "var3_none:_ChunkSizes *= *1, *1, *8, *96" &&
    ncdump -hs test_none.nc | grep -E "var2_none:_ChunkSizes *= *1, *8, *96"
  '
fi

fi
test_done
