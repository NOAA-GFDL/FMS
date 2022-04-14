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
# execute tests in the test_fms/time_manager directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

# Copy file for test.
cat <<_EOF > input_base.nml
&test_nml
test1 =.false.
test2 =.false.
test3 =.false.
test4 =.false.
test5 =.false.
test6 =.false.
test7 =.false.
test8 =.false.
test9 =.false.
test10=.false.
test11=.false.
test12=.false.
test13=.false.
test14=.false.
test15=.false.
test16=.false.
test17=.false.
test18=.false.
test19=.false.
/
_EOF

# selects next test to run through nml
testNum=0
test_next()
{
  testNum=$((testNum + 1))
  sed "s/test$testNum *=.false./test$testNum =.true./" input_base.nml > input.nml
  # test #8 must set calendar type for #9 to pass
  test $testNum -eq 9 && sed -i "s/test8 =.false./test8 =.true./" input.nml
  test_expect_success "$1" '
      mpirun -n 1 ./test_time_manager
  '
}

# run tests
test_next "set_time_i and get_time without ticks"
test_next "set_time_i and get_time with ticks"
test_next "time operators"
test_next "set_time_c"
test_next "set_date_i"
test_next "set_date_c"
test_next "increment/decrement date"
test_next "leap day cases"
test_next "days_in_month"
test_next "get_time error flag"
test_next "increment/decrement time"
test_next "negative increments"
test_next "trap for negative time"
test_next "negative seconds/ticks"
test_next "day numbering between calenders"
test_next "invalid dates"
test_next "gregorian calendar"
test_next "length_of_year"
test_next "real_to_time_type"

test_done
