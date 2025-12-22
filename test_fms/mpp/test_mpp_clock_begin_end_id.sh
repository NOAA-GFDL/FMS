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

# Chris Dupuis, Colin Gladue 06/29/20

# Set common test settings.
. ../test-lib.sh

touch input.nml
touch clock.nml
echo "&test_mpp_clock_begin_end_id_nml" > clock.nml
echo "test_number = 0" >> clock.nml
echo "/" >> clock.nml

sed -i "s/test_number = [0-9]*/test_number = 1/" clock.nml
test_expect_success "test 1" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 2/" clock.nml
test_expect_success "test 2" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 3/" clock.nml
test_expect_success "test 3" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 4/" clock.nml
test_expect_success "test 4" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 5/" clock.nml
test_expect_failure "test 5" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 6/" clock.nml
test_expect_success "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 7/" clock.nml
test_expect_success "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 8/" clock.nml
test_expect_failure "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 9/" clock.nml
test_expect_failure "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 10/" clock.nml
test_expect_failure "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 11/" clock.nml
test_expect_success "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 12/" clock.nml
test_expect_success "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 13/" clock.nml
test_expect_failure "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 14/" clock.nml
test_expect_success "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 15/" clock.nml
test_expect_failure "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 16/" clock.nml
test_expect_failure "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'

sed -i "s/test_number = [0-9]*/test_number = 17/" clock.nml
test_expect_failure "" '
    mpirun -n 1 ./test_mpp_clock_begin_end_id
'
test_done
