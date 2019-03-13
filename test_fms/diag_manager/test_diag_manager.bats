setup() {
   # Clear out any previous test files
   rm -f input.nml > /dev/null 2>&1
   rm -f diag_table > /dev/null 2>&1

   # Setup the run directory
   tnum=$( printf "%2.2d" ${BATS_TEST_NUMBER} )
   rm -f diag_test_${tnum}* > /dev/null 2>&1
   sed "s/<test_num>/${tnum}/" input.nml_base > input.nml
   ln -s diagTables/diag_table_${tnum} diag_table
}

@test "1 window, no halos" {
   run mpirun -n 1 ./test_diag_manager 
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test1.1 successful:" ]]
   [[ "$output" =~ "test1.2 successful" ]]
}

@test "multiple windows, no halos" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test2.1 successful." ]]
   [[ "$output" =~ "test2.2 successful:" ]]
}

@test "multiple windows, no halos, multiple send_data" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test3.1 successful." ]]
   [[ "$output" =~ "test3.2 successful:" ]]
}

@test "1 window, no halos, 2 time steps" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test4.1 successful." ]]
   [[ "$output" =~ "test4.2 successful:" ]]
}

@test "multiple windows, no halos, 2 time steps" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test5.1 successful." ]]
   [[ "$output" =~ "test5.2 successful:" ]]
}

@test "multiple windows, no halos, multiple send_data, 2 time steps" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test6.1 successful." ]]
   [[ "$output" =~ "test6.2 successful:" ]]
}

@test "1 window with halos" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test7.1 successful." ]]
   [[ "$output" =~ "test7.2 successful:" ]]
}

@test "1 window with halos, 2 time steps" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test8.1 successful." ]]
   [[ "$output" =~ "test8.2 successful:" ]]
}

@test "1 window, no halos, static, 1 dimension, global data, catch data too small" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test9.1 successful." ]]
   [[ "$output" =~ "test9.2 successful:" ]]
}

@test "1 window, no halos, static, 1 dimension, global data, catch data too large" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test10.1 successful." ]]
   [[ "$output" =~ "test10.2 successful:" ]]
}

@test "catch je_in missing" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test11.1 successful." ]]
   [[ "$output" =~ "test11.2 successful." ]]
}

@test "Catch duplicate field in diag_table" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -ne 0 ]
   [[ "$output" =~ "test12 successful:" ]]
}

@test "Output interval greater than runlength" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "NOTE: Potential error in diag_manager_end: dat2 NOT available, check if output interval > runlength. Netcdf fill_values are written" ]]
}

@test "Catch invalid date in register_diag_field call" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test14 successful." ]]
}

@test "OpenMP thread test" {
   sed -i -e 's/\(numthreads *= *\)[0-9]*/\12/' input.nml
   skip "TODO: Test case needs fixed"
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
}

@test "filename appendix added" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [ -e diag_test_16.g01.nc ]
}

@test "Root-mean-square" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   run nccmp -d test_data/diag_test_17.nc diag_test_17.nc
   [ "$status" -eq 0 ]
}

@test "Added attributes, and cell_measures" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   run nccmp -m test_data/diag_test_18_file1.nc diag_test_18_file1.nc
   [ "$status" -eq 0 ]
}

@test "Area and Volume same field" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -ne 0 ]
   [[ "$output" =~ "module/output_field test_mod/dat2h AREA and VOLUME CANNOT be the same variable" ]]
}

@test "Get diag_field_id, ID found and not found" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test20.1 Passes, id_dat2.EQ.id_dat2_got" ]]
   [[ "$output" =~ "test20.2 Passes, id_none_got.EQ.DIAG_FIELD_NOT_FOUND" ]]
}

@test "Add axis attributes" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   run nccmp -m test_data/diag_test_21_file1.nc diag_test_21_file1.nc
   [ "$status" -eq 0 ]
}

@test "Get 'nv' axis id" {
   run mpirun -n 1 ./test_diag_manager
   [ "$status" -eq 0 ]
   [[ "$output" =~ "test22.1 Passes: id_nv has a positive value" ]]
   [[ "$output" =~ "test22.2 Passes: Can call diag_axis_init on \"nv\" and get same ID" ]]
   run nccmp -m test_data/diag_test_22_file1.nc diag_test_22_file1.nc
   [ "$status" -eq 0 ]
   run nccmp -d test_data/diag_test_22_file1.nc diag_test_22_file1.nc
   [ "$status" -eq 0 ]
}
