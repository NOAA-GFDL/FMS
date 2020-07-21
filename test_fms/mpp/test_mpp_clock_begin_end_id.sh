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

# Chris Dupuis, Colin Gladue 06/29/20

# Set common test settings.
. ../test_common.sh


touch clock.nml
echo "&test_mpp_clock_begin_end_id_nml" > clock.nml
echo "test_number = 0" >> clock.nml
echo "/" >> clock.nml

echo
sed -i "s/test_number = [0-9]*/test_number = 1/" clock.nml
echo "Running test 1..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 1 has passed"

echo
sed -i "s/test_number = [0-9]*/test_number = 2/" clock.nml
echo "Running test 2..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 2 has passed"

echo
sed -i "s/test_number = [0-9]*/test_number = 3/" clock.nml
echo "Running test 3..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 3 has passed"

echo
sed -i "s/test_number = [0-9]*/test_number = 4/" clock.nml
echo "Running test 4..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 4 has passed"

echo
sed -i "s/test_number = [0-9]*/test_number = 5/" clock.nml
echo "Running test 5..."
err=0
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 5 was unsuccessful"
  exit 5
else
  echo "Test 5 has passed"
fi

echo
sed -i "s/test_number = [0-9]*/test_number = 6/" clock.nml
echo "Running test 6..."
err=1
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 6 was unsuccessful"
  exit 6
else
  echo "Test 6 has passed"
fi

echo
sed -i "s/test_number = [0-9]*/test_number = 7/" clock.nml
echo "Running test 7..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 7 has passed"

#echo
#sed -i "s/test_number = [0-9]*/test_number = 8/" clock.nml
#echo "Running test 8..."
#err=0
#run_test test_mpp_clock_begin_end_id 1 || err=1
#if [ "$err" -ne 1 ]; then
#  echo "ERROR: Test 8 was unsuccessful"
#  exit 8
#else
#  echo "Test 8 has passed"
#fi

echo
sed -i "s/test_number = [0-9]*/test_number = 9/" clock.nml
echo "Running test 9..."
err=0
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 9 was unsuccessful"
  exit 9
else
  echo "Test 9 has passed"
fi

echo
sed -i "s/test_number = [0-9]*/test_number = 10/" clock.nml
echo "Running test 10..."
err=0
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 10 was unsuccessful"
  exit 10
else
  echo "Test 10 has passed"
fi

echo
sed -i "s/test_number = [0-9]*/test_number = 11/" clock.nml
echo "Running test 11..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 11 has passed"

echo
sed -i "s/test_number = [0-9]*/test_number = 12/" clock.nml
echo "Running test 12..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 12 has passed"

echo
sed -i "s/test_number = [0-9]*/test_number = 13/" clock.nml
echo "Running test 13..."
err=0
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 13 was unsuccessful"
  exit 13
else
  echo "Test 13 has passed"
fi

echo
sed -i "s/test_number = [0-9]*/test_number = 14/" clock.nml
echo "Running test 14..."
run_test test_mpp_clock_begin_end_id 1
echo "Test 14 has passed"

echo
sed -i "s/test_number = [0-9]*/test_number = 15/" clock.nml
echo "Running test 15..."
err=0
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 15 was unsuccessful"
  exit 15
else
  echo "Test 15 has passed"
fi

echo
sed -i "s/test_number = [0-9]*/test_number = 16/" clock.nml
echo "Running test 16..."
err=0
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 16 was unsuccessful"
  exit 16
else
  echo "Test 16 has passed"
fi

echo
sed -i "s/test_number = [0-9]*/test_number = 17/" clock.nml
echo "Running test 17..."
err=0
run_test test_mpp_clock_begin_end_id 1 || err=1
if [ "$err" -ne 1 ]; then
  echo "ERROR: Test 17 was unsuccessful"
  exit 17
else
  echo "Test 17 has passed"
fi
