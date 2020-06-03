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

# Colin Gladue 5/28/20

# Set common test settings.
. ../test_common.sh

touch test_numb_base_ascii.nml
echo "&test_read_ascii_file_nml" > test_numb_base_ascii.nml
echo "test_numb = 0" >> test_numb_base_ascii.nml
echo "/" >> test_numb_base_ascii.nml

# Test 1
# Normal Usage
sed "s/test_numb = [0-9]/test_numb = 1/" test_numb_base_ascii.nml>test_numb_ascii.nml
# This is the ascii file we will read
cp $top_srcdir/test_fms/mpp/input_base.nml input.nml
echo "Running test 1..."
run_test test_read_ascii_file 1
echo "Test 1 has passed"

# Test 2
# get_ascii_file_num_lines not called before, fatal error
sed "s/test_numb = [0-9]/test_numb = 2/" test_numb_base_ascii.nml>test_numb_ascii.nml
echo "Running test 2..."
err=0
run_test test_read_ascii_file 1 || err=1
if [ "$err" -ne 1]; then
  echo "ERROR: Test 2 was unsuccessful"
  exit 2
else
  echo "Test 2 has passed"
fi

# Test 3
# File does not exist, fatal error
sed "s/test_numb = [0-9]/test_numb = 3/" test_numb_base_ascii.nml>test_numb_ascii.nml
echo "Running test 3..."
err=0
run_test test_read_ascii_file 1 || err=1
if [ "$err" -ne 1]; then
  echo "ERROR: Test 3 was unsuccessful"
  exit 3
else
  echo "Test 3 has passed"
fi

# Test 4
# Number of line in file is greater than size(Content(:)), fatal error
sed "s/test_numb = [0-9]/test_numb = 4/" test_numb_base_ascii.nml>test_numb_ascii.nml
touch empty.nml
echo "Running test 4..."
err=0
run_test test_read_ascii_file 1 || err=1
if [ "$err" -ne 1]; then
  echo "ERROR: Test 4 was unsuccessful"
  exit 4
else
  echo "Test 4 has passed"
fi

# Test 5
# Length of output string is too small, fatal error
sed "s/test_numb = [0-9]/test_numb = 5/" test_numb_base_ascii.nml>test_numb_ascii.nml
echo "Running test 5..."
err=0
run_test test_read_ascii_file 1 || err=1
if [ "$err" -ne 1]; then
  echo "ERROR: Test 5 was unsuccessful"
  exit 5
else
  echo "Test 5 has passed"
fi

# Test 6
# Number of lines in file does not equal to size(Content(:)), fatal error
sed "s/test_numb = [0-9]/test_numb = 6/" test_numb_base_ascii.nml>test_numb_ascii.nml
echo "Running test 6..."
err=0
run_test test_read_ascii_file 1 || err=1
if [ "$err" -ne 1]; then
  echo "ERROR: Test 6 was unsuccessful"
  exit 6
else
  echo "Test 6 has passed"
fi

# Test 7
# Normal usage, with optional PELIST argument passed in
sed "s/test_numb = [0-9]/test_numb = 7/" test_numb_base_ascii.nml>test_numb_ascii.nml
echo "Running test 7..."
run_test test_read_ascii_file 1
echo "Test 7 has passed"

# Test 8
# Normal usage, with an empty file passed in
sed "s/test_numb = [0-9]/test_numb = 8/" test_numb_base_ascii.nml>test_numb_ascii.nml
# This is the ascii file we will read
touch empty.nml
echo "Running test 8..."
run_test test_read_ascii_file 1
echo "Test 8 has passed"
