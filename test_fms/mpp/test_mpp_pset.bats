#setup() {
   ## **************
   ## Commenting out for now, until test run reliably
   ## in parallel
   ## **************
   # Clear out any previous test files
   ##rm -f input.nml > /dev/null 2>&1

   # Setup the run directory
   ##tnum=$( printf "%2.2d" ${BATS_TEST_NUMBER} )
   ##rm -f diag_test_${tnum}* > /dev/null 2>&1
   ##sed "s/<test_num>/${tnum}/" input.nml_base > input.nml
#}

teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "Test 1: Tests how to distribute allocatable arrays" {
    skip "Test unreliable for parallel tests"
    run mpirun -n 2 ./test_mpp_pset
    [ "$status" -eq 0 ]
}

@test "Test 2: Tests how to distribute automatic arrays" {
    skip "Test unreliable for parallel tests"
    run mpirun -n 2 ./test_mpp_pset
    [ "$status" -eq 0 ]
}
