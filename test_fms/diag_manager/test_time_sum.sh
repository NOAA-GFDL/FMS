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

# Copyright (c) 2019-2020 Ed Hartnett, Seth Underwood

# Set common test settings.
. ../test-lib.sh

if [ -z "${skipflag}" ]; then
# create and enter directory for in/output files
output_dir

#TODO replace with yaml diag_table and set diag_manager_nml::use_modern_diag=.true.
cat <<_EOF > diag_table
test_sum
2 1 1 0 0 0

"test_sum",      6,  "hours", 1, "hours", "time"
"test_sum_regional",      6,  "hours", 1, "hours", "time"

"ocn_mod", "var0", "var0_sum", "test_sum", "all", "sum", "none", 2
"ocn_mod", "var1", "var1_sum", "test_sum", "all", "sum", "none", 2
"ocn_mod", "var2", "var2_sum", "test_sum", "all", "sum", "none", 2
"ocn_mod", "var3", "var3_sum", "test_sum", "all", "sum", "none", 2

"ocn_mod", "var3", "var3_Z", "test_sum", "all", "sum", "-1 -1 -1 -1 2. 3.", 2

"ocn_mod", "var3", "var3_sum", "test_sum_regional", "all", "sum", "78. 81. 78. 81. 2. 3.", 2 #chosen by MKL
_EOF

my_test_count=1
printf "&test_reduction_methods_nml \n test_case = 0 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "sum" reduction method (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'

my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 0 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_failure "Checking answers for the "sum" reduction method, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'

my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 0 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_failure "Checking answers for the "sum" reduction method, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'

export OMP_NUM_THREADS=2
my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method with openmp (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_failure "Checking answers for the "sum" reduction method with openmp (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'

my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 1 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_failure "Checking answers for the "sum" reduction method with openmp, logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'

my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 1 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_failure "Checking answers for the "sum" reduction method with openmp, real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'
export OMP_NUM_THREADS=1

my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method with halo output (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_success "Checking answers for the "sum" reduction method with halo output (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'

my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 2 \n mask_case = 1 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_failure "Checking answers for the "sum" reduction method with halo output with logical mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'

my_test_count=`expr $my_test_count + 1`
printf "&test_reduction_methods_nml \n test_case = 2 \n mask_case = 2 \n \n/" | cat > input.nml
test_expect_success "Running diag_manager with "sum" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 6 ../test_reduction_methods
'
test_expect_failure "Checking answers for the "sum" reduction method with halo output with real mask (test $my_test_count)" '
  mpirun -n 1 ../check_time_sum
'
fi
test_done