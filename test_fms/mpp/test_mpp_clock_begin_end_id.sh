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
. ../test-lib.sh

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
