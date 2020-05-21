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

# Colin Gladue 04/27/2020

# Set common test settings.
. ../test_common.sh

touch test_numb_base.nml
echo "&test_read_input_nml_nml" > test_numb_base.nml
echo "test_numb = 0" >> test_numb_base.nml
echo "/" >> test_numb_base.nml

# Test 1
sed "s/test_numb = [0-9]/test_numb = 1/" test_numb_base.nml > test_numb.nml
cp $top_srcdir/test_fms/mpp/input_base.nml input.nml
run_test test_read_input_nml 1
if [ $? = 0 ]; then # Checks if running the subroutine causes an error or not
  awk '{ sub(/^[ \t]+/, ""); print }' input.nml > trimmed_input_test1.tst
  awk '{ sub(/^[ \t]+/, ""); print }' logfile.000000.out > trimmed_log_test1.tst
  sort trimmed_input_test1.tst > sorted_input_test1.tst
  sort trimmed_log_test1.tst > sorted_log_test1.tst
  input_var1=$(comm -12 sorted_input_test1.tst sorted_input_test1.tst) # Done this way to achieve same formatting as next line
  log_var1=$(comm -12 sorted_log_test1.tst sorted_log_test1.tst) # Done this way to achieve same formatting as next line
  incommon_var1=$(comm -12 sorted_input_test1.tst sorted_log_test1.tst)
  if [ "$input_var1" = "$incommon_var1" ]; then # Checks if the logfile contains all of the input nml
    err=0
    grep -n "READ_INPUT_NML: input.nml" logfile.000000.out|| err=1
    grep -n "READ_INPUT_NML: unknown" logfile.000000.out|| err=1
    if [ "$err" != 1 ]; then # Checks if the logfile lists the version and filename
      echo "Test 1 has passed"
    else
      echo "ERROR: Test 1 was unsuccessful. Version or filename not correctly written."
      exit 31
    fi
  else
    echo "ERROR: Test 1 was unsuccessful. Log did not contain input.nml"
    exit 21
  fi
else
  echo "ERROR: Test 1 was unsuccessful. The read_input_nml subroutine failed to execute"
  exit 11
fi

# Test 2
sed "s/test_numb = [0-9]/test_numb = 2/" test_numb_base.nml > test_numb.nml
sed "s/1/2/" $top_srcdir/test_fms/mpp/input_base.nml > input_alternative.nml
run_test test_read_input_nml 1
if [ $? = 0 ]; then
  awk '{ sub(/^[ \t]+/, ""); print }' input_alternative.nml > trimmed_input_test2.tst
  awk '{ sub(/^[ \t]+/, ""); print }' logfile.000000.out > trimmed_log_test2.tst
  sort trimmed_input_test2.tst > sorted_input_test2.tst
  sort trimmed_log_test2.tst > sorted_log_test2.tst
  input_var2=$(comm -12 sorted_input_test2.tst sorted_input_test2.tst) # Done this way to achieve same formatting as next line
  log_var2=$(comm -12 sorted_log_test2.tst sorted_log_test2.tst) # Done this way to achieve same formatting as next line
  incommon_var2=$(comm -12 sorted_input_test2.tst sorted_log_test2.tst)
  if [ "$input_var2" = "$incommon_var2" ]; then
    err=0
    grep -n "READ_INPUT_NML: input_alternative.nml" logfile.000000.out|| err=1
    grep -n "READ_INPUT_NML: unknown" logfile.000000.out|| err=1
    if [ "$err" != 1 ]; then # Checks if the logfile lists the version and filename
      echo "Test 2 has passed"
    else
      echo "ERROR: Test 2 was unsuccessful. Version or filename not correctly written."
      exit 32
    fi
  else
    echo "ERROR: Test 2 was unsuccessful. Log did not contain input_alternative.nml"
    exit 22
  fi
else
  echo "ERROR: Test 2 was unsuccessful. The read_input_nml subroutine failed to execute"
  exit 12
fi

# Test 3
sed "s/test_numb = [0-9]/test_numb = 3/" test_numb_base.nml > test_numb.nml
err=0
run_test test_read_input_nml 1 || err=1
if [ "$err" != 1 ]; then
  echo "ERROR: Test 3 was unsuccessful."
  exit 13
else
   echo "Test 3 has passed"
fi

# Test 4
sed "s/test_numb = [0-9]/test_numb = 4/" test_numb_base.nml > test_numb.nml
rm input.nml
touch input.nml # Achieve a blank namelist to be read
run_test test_read_input_nml 1
if [ $? = 0 ]; then
  err=0
  grep -n "READ_INPUT_NML: input.nml" logfile.000000.out|| err=1
  grep -n "READ_INPUT_NML: unknown" logfile.000000.out|| err=1
  if [ "$err" != 1 ]; then # Checks if the logfile lists the version and filename
    echo "Test 4 has passed"
  else
    echo "ERROR: Test 4 was unsuccessful. Version or filename not correctly written."
    exit 34
  fi
else
  echo "ERROR: Test 4 was unsuccessful. Subroutine failed to execute"
  exit 14
fi
