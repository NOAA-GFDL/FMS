#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/data_override directory.

# Ed Hartnett 11/26/19

# Set common test settings.
. ../test_common.sh

copy_files()
{

    tnum=$( printf "%2.2d" ${1} )
    rm -f diag_test_${tnum}* > /dev/null 2>&1
    sed "s/<test_num>/${tnum}/" $top_srcdir/test_fms/diag_manager/input.nml_base > input.nml
    ln -f -s $top_srcdir/test_fms/diag_manager/diagTables/diag_table_${tnum} diag_table
}

run_test()
{
    echo ${2}
    copy_files ${1}
    if test "x${3}" = "xfail"; then
        if mpirun -n 1 ./test_diag_manager; then exit 1; fi
    else
        mpirun -n 1 ./test_diag_manager
    fi
}

rm -f input.nml diag_table
run_test 1 "Test 1: Data array is too large in x and y direction"
run_test 2 "Test 2: Data array is too large in x direction"
run_test 3 "Test 3: Data array is too large in y direction"
run_test 4 "Test 4: Data array is too small in x and y direction, checks for 2 time steps"
run_test 5 "Test 5: Data array is too small in x directions, checks for 2 time steps"
run_test 6 "Test 6: Data array is too small in y direction, checks for 2 time steps"
run_test 7 "Test 7: Data array is too large in x and y, with halos, 2 time steps"
run_test 8 "Test 8: Data array is too small in x and y, with halos, 2 time steps"
run_test 9 "Test 9: Data array is too small, 1D, static global data"
run_test 10 "Test 10: Data array is too large, 1D, static global data"
run_test 11 "Test 11: Missing je_in as an input"
run_test 12 "Test 12: Catch duplicate field in diag_table" "fail"
run_test 13 "Test 13: Output interval greater than runlength"
run_test 14 "Test 14: Catch invalid date in register_diag_field call"
run_test 15 "Test 15: OpenMP thread test"
run_test 16 "Test 16: Filename appendix added"
run_test 17 "Test 17: Root-mean-square"
run_test 18 "Test 18: Added attributes, and cell_measures"
run_test 19 "Test 19: Area and Volume same field" "fail"
run_test 20 "Test 20: Get diag_field_id, ID found and not found"
run_test 21 "Test 21: Add axis attributes"
run_test 22 "Test 22: Get 'nv' axis id"
run_test 23 "Test 23: Unstructured grid"
