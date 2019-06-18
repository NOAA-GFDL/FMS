setup() {
   # Clear out any previous test files
   rm -f input.nml > /dev/null 2>&1

   # Setup the run directory
   tnum=$( printf "%2.2d" ${BATS_TEST_NUMBER} )
   sed "s/<test_num>/${tnum}/" input.nml_base > input.nml
}

@test "Test 1: Cubic-Grid" {
   run mpirun -n 6 ./test_data_override
   [ "$status" -eq 0 ]
}

@test "Test 2: Latlon-Grid" {
   run mpirun -n 6 ./test_data_override
   [ "$status" -eq 0 ]
}

