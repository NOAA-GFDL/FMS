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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Eric Stofferahn 07/15/2020

# Set common test settings.
. ../test_common.sh

skip_test="no"

# Copy file for test.
cp $top_srcdir/test_fms/mpp/base_ascii_5 ascii_5
cp $top_srcdir/test_fms/mpp/base_ascii_25 ascii_25
cp $top_srcdir/test_fms/mpp/base_ascii_0 ascii_0
cp $top_srcdir/test_fms/mpp/base_ascii_skip ascii_skip
cp $top_srcdir/test_fms/mpp/base_ascii_long ascii_long

# Set up namelist to carry test_number.
touch test_numb_base2.nml
echo "&test_mpp_get_ascii_lines_nml" > test_numb_base2.nml
echo "test_number = 0" >> test_numb_base2.nml
echo "/" >> test_numb_base2.nml

for tst in 1 2 3 4
do
  sed "s/test_number = [0-9]/test_number = ${tst}/" test_numb_base2.nml > test_numb2.nml
  echo "Running test ${tst}..."
  run_test test_mpp_get_ascii_lines 2 $skip_test
  echo "Test ${tst} has passed"
done

sed "s/test_number = [0-9]/test_number = 5/" test_numb_base2.nml > test_numb2.nml
echo "Running test 5..."
run_test test_mpp_get_ascii_lines 2 $skip_test || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 5 was unsuccessful."
  exit 5
else
   echo "Test 5 has passed"
fi
