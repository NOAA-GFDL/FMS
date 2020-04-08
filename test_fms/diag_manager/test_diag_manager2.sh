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
# execute tests in the test_fms/data_override directory.

# Ed Hartnett 11/26/19

# Set common test settings.
. ../test_common.sh

setup_test()
# Function sets up and runs the test
# Inputs:
# {1} Test number
# {2} Description of test
# {3} Set to "skip"

{
    echo ${2}
    tnum=$( printf "%2.2d" ${1} )
    rm -f diag_test_${tnum}* > /dev/null 2>&1
    sed "s/<test_num>/${tnum}/" $top_srcdir/test_fms/diag_manager/input.nml_base > input.nml
    ln -f -s $top_srcdir/test_fms/diag_manager/diagTables/diag_table_${tnum} diag_table

    run_test test_diag_manager 1 ${3}

}

rm -f input.nml diag_table
setup_test 1 "Test 1: Data array is too large in x and y direction"
setup_test 2 "Test 2: Data array is too large in x direction"
setup_test 3 "Test 3: Data array is too large in y direction"
setup_test 4 "Test 4: Data array is too small in x and y direction, checks for 2 time steps"
setup_test 5 "Test 5: Data array is too small in x directions, checks for 2 time steps"
setup_test 6 "Test 6: Data array is too small in y direction, checks for 2 time steps"
setup_test 7 "Test 7: Data array is too large in x and y, with halos, 2 time steps"
setup_test 8 "Test 8: Data array is too small in x and y, with halos, 2 time steps"
setup_test 9 "Test 9: Data array is too small, 1D, static global data"
setup_test 10 "Test 10: Data array is too large, 1D, static global data"
setup_test 11 "Test 11: Missing je_in as an input"
setup_test 12 "Test 12: Catch duplicate field in diag_table"
setup_test 13 "Test 13: Output interval greater than runlength"
setup_test 14 "Test 14: Catch invalid date in register_diag_field call"
setup_test 15 "Test 15: OpenMP thread test"
setup_test 16 "Test 16: Filename appendix added"
setup_test 17 "Test 17: Root-mean-square"
setup_test 18 "Test 18: Added attributes, and cell_measures"
setup_test 19 "Test 19: Area and Volume same field"
setup_test 20 "Test 20: Get diag_field_id, ID found and not found"
setup_test 21 "Test 21: Add axis attributes"
setup_test 22 "Test 22: Get 'nv' axis id"
setup_test 23 "Test 23: Unstructured grid"
