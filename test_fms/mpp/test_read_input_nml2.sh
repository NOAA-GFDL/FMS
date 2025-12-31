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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Colin Gladue 05/22/2020
# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

# create and enter directory for in/output
output_dir

touch test_numb_base.nml
echo "&test_read_input_nml_nml" > test_numb_base.nml
echo "test_numb = 0" >> test_numb_base.nml
echo "/" >> test_numb_base.nml
cp test_numb_base.nml input.nml

# Test 1
sed "s/test_numb = [0-9]/test_numb = 1/" test_numb_base.nml > test_numb.nml
test_expect_success "read input nml" '
    mpirun -n 1 ../test_read_input_nml
'

# Test 2
sed "s/test_numb = [0-9]/test_numb = 2/" test_numb_base.nml > test_numb.nml
cp input.nml input_alternative.nml
test_expect_success "read input nml with file name" '
    mpirun -n 1 ../test_read_input_nml
'

# Test 3
sed "s/test_numb = [0-9]/test_numb = 3/" test_numb_base.nml > test_numb.nml
test_expect_failure "failure caught on invalid nml" '
    mpirun -n 1 ../test_read_input_nml
'

# Test 4
sed "s/test_numb = [0-9]/test_numb = 4/" test_numb_base.nml > test_numb.nml
touch input_blank.nml # blank namelist to be read
test_expect_success "read empty nml" '
    mpirun -n 1 ../test_read_input_nml
'

test_done
