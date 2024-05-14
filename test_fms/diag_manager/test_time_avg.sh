#!/bin/sh

#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# Set common test settings.
. ../test-lib.sh

if [ -z "${skipflag}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.yaml
title: test_avg
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_avg
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    output_name: var0_avg
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var1
    output_name: var1_avg
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var2
    output_name: var2_avg
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3_avg
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var4
    output_name: var4_avg
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3_Z
    reduction: average
    zbounds: 2. 3.
    kind: r4
  - module: ocn_mod
    var_name: IOnASphere
    reduction: average
    kind: r4
- file_name: test_avg_regional
  time_units: hours
  unlimdim: time
  sub_region:
  - grid_type: latlon
    corner1: 78. 78.
    corner2: 78. 78.
    corner3: 81. 81.
    corner4: 81. 81.
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var3
    output_name: var3_avg
    reduction: average
    zbounds: 2. 3.
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc

# tests with no mask, no openmp
my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 0 \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n &test_reduction_methods_nml \n test_case = 0 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n &test_reduction_methods_nml \n test_case = 0 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

# openmp tests

export OMP_NUM_THREADS=2
my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method with openmp (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method with openmp (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

# halo output and mask tests

export OMP_NUM_THREADS=1

# This is the corner case where the number of openmp threads is 1 but the number of
# atmosphere blocks is not set 1!
my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "none" reduction method with blocking but no threads (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "none" reduction method with blocking but no threads (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method with halo output (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method with halo output (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "avg" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "avg" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_avg
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n mix_snapshot_average_fields = .true. \n /" | cat > input.nml
test_expect_failure "Running diag_manager with with mix_snapshot_average_fields = .true. (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'

cat <<_EOF > diag_table.yaml
title: test_avg
base_date: 2 1 1 0 0 0
diag_files:
- file_name:  test_avg
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    output_name: var0_avg
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var0
    output_name: var0_none
    reduction: none
    kind: r4
_EOF

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n /" | cat > input.nml
test_expect_failure "Running diag_manager with with a file with instantaneous and averaged output (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'

cat <<_EOF > diag_table.yaml
title: test_avg
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_failure
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var2missing
    reduction: average
    kind: r4
_EOF

  my_test_count=`expr $my_test_count + 1`
  test_expect_failure "Fail if passing in missing_values without masking them (test $my_test_count)" '
    mpirun -n 6 ../test_reduction_methods
  '
fi
test_done
