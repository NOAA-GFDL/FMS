#!/bin/sh

#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL Flexible Modeling System (FMS).
#
# FMS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FMS is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Colin Gladue 05/22/2020
# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

# create and enter directory for in/output
output_dir

touch input.nml
touch test_numb_base.nml
echo "&test_read_input_nml_nml" > test_numb_base.nml
echo "test_numb = 0" >> test_numb_base.nml
echo "/" >> test_numb_base.nml

# Test 1
sed "s/test_numb = [0-9]/test_numb = 1/" test_numb_base.nml > test_numb.nml
test_expect_success "read input nml" '
    mpirun -n 1 ../test_read_input_nml
'

# Test 2
sed "s/test_numb = [0-9]/test_numb = 2/" test_numb_base.nml > test_numb.nml
sed "s/1/2/" $top_srcdir/test_fms/mpp/input_base.nml > input_alternative.nml
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
