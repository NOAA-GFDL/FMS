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

SKIP_TESTS="test_time_diurnal.[4-20]"

if [ -z "${skipflag}" ]; then
# create and enter directory for in/output files
output_dir

# diurnal specific test
# aims to reproduce files from diurnal test (test_diag_manager_time)

cat <<_EOF > diag_table.yaml
title: test_diurnal
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_diurnal
  time_units: hours 
  unlimdim: time
  freq: 1 months 
  varlist:
  - module: ocn_mod 
    var_name: sst
    output_name: sst
    reduction: diurnal3
    kind: r4
  - module: ocn_mod 
    var_name: ice 
    output_name: ice 
    reduction: diurnal3
    kind: r4
_EOF

printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n" > input.nml

test_expect_success "simple hourly diurnal test" '
  mpirun -n 1 ../test_diag_diurnal
'

test_expect_success "checking results for simple hourly diurnal test" '
  mpirun -n 1 ../check_time_diurnal
' 

# end here for now
test_done
exit

# this passes but doesn't seem to be using the new code as it should
test_expect_success "old diurnal test with modern_diag" '
  mpirun -n 1 ../test_diag_manager_time
'


# use old code for now
#rm input.nml && touch input.nml
#cat <<_EOF > diag_table
#test_diag_manager
#2 1 1 0 0 0
#
##output files
#"test_diurnal",         1, "hours",   1, "hours", "time"
#
##output variables
# "test_diag_manager_mod", "sst", "sst", "test_diurnal",  "all", "diurnal4", "none", 2
# "test_diag_manager_mod", "ice", "ice", "test_diurnal",  "all", "diurnal4", "none", 2
#_EOF
#test_expect_success "old diurnal test (with new code)" '
#  mpirun -n 1 ../test_diag_manager_time
#'

## just run the very basic test for now
test_done
exit

# same as the other reduction tests

cat <<_EOF > diag_table.yaml
title: test_diurnal
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_diurnal
  time_units: hours
  unlimdim: time
  freq: 6 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    output_name: var0
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var1
    output_name: var1
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var2
    output_name: var2
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var4
    output_name: var4
    reduction: diurnal3
    kind: r4
  - module: ocn_mod
    var_name: var3
    output_name: var3_Z
    reduction: diurnal3
    zbounds: 2. 3.
    kind: r4
- file_name: test_regional
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
    output_name: var3
    reduction: diurnal3
    zbounds: 2. 3.
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc

# tests with no mask, no openmp
my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 0 \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'


my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n &test_reduction_methods_nml \n test_case = 0 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n &test_reduction_methods_nml \n test_case = 0 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'

# openmp tests

export OMP_NUM_THREADS=2
my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method with openmp (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method with openmp (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 1 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'

# halo output and mask tests

export OMP_NUM_THREADS=1

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method with halo output (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method with halo output (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'

my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n / \n&test_reduction_methods_nml \n test_case = 2 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "diurnal" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "diurnal" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_diurnal
'
fi
test_done
