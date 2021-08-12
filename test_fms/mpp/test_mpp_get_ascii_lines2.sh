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
# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

# create and enter directory for in/output
output_dir

# create ascii files
cat <<_EOF > ascii_5
"this is an ascii file with 5 lines"
"it will contain commas inside quotes", "as well as outside quotes"
"it will not have the same number of fields on every line"
some lines will not have quotes
"there might be a line with the character string \n"
_EOF

cat <<_EOF > ascii_long
"this is an ascii file with 5 lines"
""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes""it will contain commas inside quotes", "as well as outside quotes"it will contain commas inside quotes", "as well as outside quotes"
"it will not have the same number of fields on every line"
some lines will not have quotes
"there might be a line with the character string \n"
_EOF

cat <<_EOF > ascii_skip
"this is an ascii file with 5 lines"
"it will contain commas inside quotes", "as well as outside quotes"
"it will not have the same number of fields on every line"

"there might be a blank line beforehand"
_EOF

touch ascii_25
for i in 1 2 3 4 5
do
cat ascii_5 >> ascii_25
done

echo "" > ascii_0

# Set up namelist to carry test_number.
touch input.nml
touch test_numb_base2.nml
echo "&test_mpp_get_ascii_lines_nml" > test_numb_base2.nml
echo "test_number = 0" >> test_numb_base2.nml
echo "/" >> test_numb_base2.nml

# run tests
sed "s/test_number = [0-9]/test_number = 1/" test_numb_base2.nml > test_numb2.nml
test_expect_success "5 lines" '
    mpirun -n 2 ../test_mpp_get_ascii_lines
'
sed "s/test_number = [0-9]/test_number = 2/" test_numb_base2.nml > test_numb2.nml
test_expect_success "25 lines" '
    mpirun -n 2 ../test_mpp_get_ascii_lines
'
sed "s/test_number = [0-9]/test_number = 3/" test_numb_base2.nml > test_numb2.nml
test_expect_success "0 lines" '
    mpirun -n 2 ../test_mpp_get_ascii_lines
'
sed "s/test_number = [0-9]/test_number = 4/" test_numb_base2.nml > test_numb2.nml
test_expect_success "blank line" '
    mpirun -n 2 ../test_mpp_get_ascii_lines
'
sed "s/test_number = [0-9]/test_number = 5/" test_numb_base2.nml > test_numb2.nml
test_expect_failure "failure caught from long line" '
    mpirun -n 2 ../test_mpp_get_ascii_lines
'
test_done
