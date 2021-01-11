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

# Loop through tests
testNum=1
while [ $testNum -le 19 ]
do
    sed "s/test$testNum *=.false./test$testNum =.true./" input_base.nml > input.nml
    test_expect_success "time manager test #$testNum" '
        mpirun -n 1 ./test_time_manager
    '
    testNum=$((testNum + 1))
done

test_done
