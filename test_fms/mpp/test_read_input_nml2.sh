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

# Set common test settings.
. ../test_common.sh

touch test_numb_base.nml
echo "&test_read_input_nml_nml" > test_numb_base.nml
echo "test_numb = 0" >> test_numb_base.nml
echo "/" >> test_numb_base.nml

# Test 1
sed "s/test_numb = [0-9]/test_numb = 1/" test_numb_base.nml > test_numb.nml
cp $top_srcdir/test_fms/mpp/input_base.nml input.nml
echo "Running test 1..."
run_test test_read_input_nml 1
echo "Test 1 has passed"

# Test 2
sed "s/test_numb = [0-9]/test_numb = 2/" test_numb_base.nml > test_numb.nml
sed "s/1/2/" $top_srcdir/test_fms/mpp/input_base.nml > input_alternative.nml
echo "Running test 2..."
run_test test_read_input_nml 1
echo "Test 2 has passed"

# Test 3
sed "s/test_numb = [0-9]/test_numb = 3/" test_numb_base.nml > test_numb.nml
echo "Running test 3..."
run_test test_read_input_nml 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 3 was unsuccessful."
  exit 3
else
   echo "Test 3 has passed"
fi

# Test 4
sed "s/test_numb = [0-9]/test_numb = 4/" test_numb_base.nml > test_numb.nml
touch input_blank.nml # Achieve a blank namelist to be read
echo "Running test 4..."
run_test test_read_input_nml 1
echo "Test 4 has passed"
