setup() {
   # Clear out any previous test files
   rm -f input.nml > /dev/null 2>&1

   # Setup the run directory
   tnum=$( printf "%2.2d" ${BATS_TEST_NUMBER} )
   rm -f diag_test_${tnum}* > /dev/null 2>&1
   sed "s/<test_num>/${tnum}/" input.nml_base > input.nml
}

@test "Test 1: Tests how to distribute allocatable arrays" {
    run mpirun -n 6 ./test_mpp_pset
    [ "$status" -eq 0 ]
}

@test "Test 2: Tests how to distribute automatic arrays" {
    run mpirun -n 6 ./test_mpp_pset
    [ "$status" -eq 0 ]
}
