setup() {
   # Clear out any previous test files
   rm -f input.nml > /dev/null 2>&1
   rm -f diag_table > /dev/null 2>&1

   # Setup the run directory
   tnum=$( printf "%2.2d" ${BATS_TEST_NUMBER} )
   rm -f diag_test_${tnum}* > /dev/null 2>&1
<<<<<<< HEAD
   sed "s/<test_num>/${tnum}/" ${srcdir}/input.nml_base > input.nml
   ln -s ${srcdir}/diagTables/diag_table_${tnum} diag_table
}

@test "1 window, no halos" {
   skip "TODO: figure out why this fails."
   run mpirun -n 1 ./test_diag_manager 
=======
   sed "s/<test_num>/${tnum}/" $srcdir/test_fms/diag_manager/input.nml_base > input.nml
   ln -s $srcdir/test_fms/diag_manager/diagTables/diag_table_${tnum} diag_table
}

teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "Test 1: Data array is too large in x and y direction" {
   run mpirun -n 1 ./test_diag_manager
>>>>>>> master
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test1.1 successful:" ]]
   [[ "$output" =~ "test1.2 successful" ]]
}

<<<<<<< HEAD
@test "multiple windows, no halos" {
   skip "TODO: figure out why this fails."
=======

@test "Test 2: Data array is too large in x direction" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test2.1 successful." ]]
   [[ "$output" =~ "test2.2 successful:" ]]
}

<<<<<<< HEAD
@test "multiple windows, no halos, multiple send_data" {
   skip "TODO: figure out why this fails."
=======
@test "Test 3: Data array is too large in y direction" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test3.1 successful." ]]
   [[ "$output" =~ "test3.2 successful:" ]]
}

<<<<<<< HEAD
@test "1 window, no halos, 2 time steps" {
   skip "TODO: figure out why this fails."
=======
@test "Test 4: Data array is too small in x and y direction, checks for 2 time steps" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test4.1 successful." ]]
   [[ "$output" =~ "test4.2 successful:" ]]
}

<<<<<<< HEAD
@test "multiple windows, no halos, 2 time steps" {
   skip "TODO: figure out why this fails."
=======
@test "Test 5: Data array is too small in x directions, checks for 2 time steps" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test5.1 successful." ]]
   [[ "$output" =~ "test5.2 successful:" ]]
}

<<<<<<< HEAD
@test "multiple windows, no halos, multiple send_data, 2 time steps" {
   skip "TODO: figure out why this fails."
=======
@test "Test 6: Data array is too small in y direction, checks for 2 time steps" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test6.1 successful." ]]
   [[ "$output" =~ "test6.2 successful:" ]]
}

<<<<<<< HEAD
@test "1 window with halos" {
   skip "TODO: figure out why this fails."
=======
@test "Test 7: Data array is too large in x and y, with halos, 2 time steps" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test7.1 successful." ]]
   [[ "$output" =~ "test7.2 successful:" ]]
}

<<<<<<< HEAD
@test "1 window with halos, 2 time steps" {
   skip "TODO: figure out why this fails."
=======
@test "Test 8: Data array is too small in x and y, with halos, 2 time steps" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test8.1 successful." ]]
   [[ "$output" =~ "test8.2 successful:" ]]
}

<<<<<<< HEAD
@test "1 window, no halos, static, 1 dimension, global data, catch data too small" {
   skip "TODO: figure out why this fails."
=======
@test "Test 9: Data array is too small, 1D, static global data" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test9.1 successful." ]]
   [[ "$output" =~ "test9.2 successful:" ]]
}

<<<<<<< HEAD
@test "1 window, no halos, static, 1 dimension, global data, catch data too large" {
   skip "TODO: figure out why this fails."
=======
@test "Test 10: Data array is too large, 1D, static global data" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test10.1 successful." ]]
   [[ "$output" =~ "test10.2 successful:" ]]
}

<<<<<<< HEAD
@test "catch je_in missing" {
   skip "TODO: figure out why this fails."
=======
@test "Test 11: Missing je_in as an input" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test11.1 successful." ]]
   [[ "$output" =~ "test11.2 successful." ]]
}

<<<<<<< HEAD
@test "Catch duplicate field in diag_table" {
   skip "TODO: figure out why this fails."
=======
@test "Test 12: Catch duplicate field in diag_table" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -ne 0 ]
   [[ "$output" =~ "test12 successful:" ]]
}

<<<<<<< HEAD
@test "Output interval greater than runlength" {
   skip "TODO: figure out why this fails."
=======
@test "Test 13: Output interval greater than runlength" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "NOTE: Potential error in diag_manager_end: dat2 NOT available, check if output interval > runlength. Netcdf fill_values are written" ]]
}

<<<<<<< HEAD
@test "Catch invalid date in register_diag_field call" {
   skip "TODO: figure out why this fails."
=======
@test "Test 14: Catch invalid date in register_diag_field call" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test14 successful." ]]
}

@test "Test 15: OpenMP thread test" {
   skip
   sed -i -e 's/\(numthreads *= *\)[0-9]*/\12/' input.nml
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
}

<<<<<<< HEAD
@test "filename appendix added" {
   skip "TODO: figure out why this fails."
=======
@test "Test 16: Filename appendix added" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [ -e diag_test_16.g01.nc ]
}

<<<<<<< HEAD
@test "Root-mean-square" {
   skip "TODO: figure out why this fails."
=======
@test "Test 17: Root-mean-square" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
#   run ncdump diag_test_17.nc
#   echo "$output" > diag_test_17.out
#   run cmp diag_test_17.out test_data/diag_test_17.out
#   [ "$status" -eq 0 ]
}

<<<<<<< HEAD
@test "Added attributes, and cell_measures" {
   skip "TODO: figure out why this fails."
=======
@test "Test 18: Added attributes, and cell_measures" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
#   run ncdump diag_test_18_file1.nc
#   echo "$output" > diag_test_18_file1.out
#   run cmp test_data/diag_test_18_file1.out diag_test_18_file1.out
#   [ "$status" -eq 0 ]
}

<<<<<<< HEAD
@test "Area and Volume same field" {
   skip "TODO: figure out why this fails."
=======
@test "Test 19: Area and Volume same field" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -ne 0 ]
   [[ "$output" =~ "module/output_field test_mod/dat2h AREA and VOLUME CANNOT be the same variable" ]]
}

<<<<<<< HEAD
@test "Get diag_field_id, ID found and not found" {
   skip "TODO: figure out why this fails."
=======
@test "Test 20: Get diag_field_id, ID found and not found" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test20.1 Passes, id_dat2.EQ.id_dat2_got" ]]
   [[ "$output" =~ "test20.2 Passes, id_none_got.EQ.DIAG_FIELD_NOT_FOUND" ]]
}

<<<<<<< HEAD
@test "Add axis attributes" {
   skip "TODO: figure out why this fails."
=======
@test "Test 21: Add axis attributes" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
#   run ncdump diag_test_21_file1.nc
#   echo "$output" > diag_test_21_file1.out
#   run cmp test_data/diag_test_21_file1.out diag_test_21_file1.out
#   [ "$status" -eq 0 ]
}

<<<<<<< HEAD
@test "Get 'nv' axis id" {
   skip "TODO: figure out why this fails."
=======
@test "Test 22: Get 'nv' axis id" {
>>>>>>> master
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test22.1 Passes: id_nv has a positive value" ]]
   [[ "$output" =~ "test22.2 Passes: Can call diag_axis_init on \"nv\" and get same ID" ]]
#   run ncdump diag_test_22_file1.nc
#   echo "$output" > diag_test_22_file1.out
#   run cmp test_data/diag_test_22_file1.out diag_test_22_file1.out
#   [ "$status" -eq 0 ]
}

@test "Test 23: Unstructured grid" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
}
