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
# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

# create and enter directory for in/output
output_dir

# create input namelist
cat <<_EOF > input.nml
&test_mpp_get_ascii_lines_nml
test_number = <test_num>
/
_EOF

touch test_numb_base_ascii.nml
echo "&test_read_ascii_file_nml" > test_numb_base_ascii.nml
echo "test_numb = 0" >> test_numb_base_ascii.nml
echo "/" >> test_numb_base_ascii.nml

# Test 1
# Normal Usage
sed "s/test_numb = [0-9]/test_numb = 1/" test_numb_base_ascii.nml>test_numb_ascii.nml
test_expect_success "normal ascii usage" '
    mpirun -n 1 ../test_read_ascii_file
'

# Test 2
# get_ascii_file_num_lines not called before, fatal error
sed "s/test_numb = [0-9]/test_numb = 2/" test_numb_base_ascii.nml>test_numb_ascii.nml
test_expect_failure "failure caught if get_ascii_file_num_lines not called before" '
    mpirun -n 1 ../test_read_ascii_file
'

# Test 3
# File does not exist, fatal error
sed "s/test_numb = [0-9]/test_numb = 3/" test_numb_base_ascii.nml>test_numb_ascii.nml
test_expect_failure "failure caught if file does not exist" '
    mpirun -n 1 ../test_read_ascii_file
'

# Test 4
# Number of line in file is greater than size(Content(:)), fatal error
sed "s/test_numb = [0-9]/test_numb = 4/" test_numb_base_ascii.nml>test_numb_ascii.nml
echo "" > empty.nml
test_expect_failure "failure caught from too few input lines" '
    mpirun -n 1 ../test_read_ascii_file
'
# Test 5
# Length of output string is too small, fatal error
sed "s/test_numb = [0-9]/test_numb = 5/" test_numb_base_ascii.nml>test_numb_ascii.nml
test_expect_failure "failure caught from too small output string" '
    mpirun -n 1 ../test_read_ascii_file
'

# Test 6
# Number of lines in file does not equal to size(Content(:)), fatal error
sed "s/test_numb = [0-9]/test_numb = 6/" test_numb_base_ascii.nml>test_numb_ascii.nml
test_expect_failure "failure caught from mismatching numbers of lines" '
    mpirun -n 1 ../test_read_ascii_file
'

# Test 7
# Normal usage, with optional PELIST argument passed in
sed "s/test_numb = [0-9]/test_numb = 7/" test_numb_base_ascii.nml>test_numb_ascii.nml
test_expect_success "normal ascii usage with PELIST" '
    mpirun -n 1 ../test_read_ascii_file
'

# Test 8
# Normal usage, with an empty file passed in
sed "s/test_numb = [0-9]/test_numb = 8/" test_numb_base_ascii.nml>test_numb_ascii.nml
touch empty.nml
test_expect_success "normal ascii usage with empty file" '
    mpirun -n 1 ../test_read_ascii_file
'

test_done
