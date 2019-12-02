setup() {
   # Clear out any previous test files
   rm -f input.nml > /dev/null 2>&1

   # Setup the run directory
   tnum=$( printf "%2.2d" ${BATS_TEST_NUMBER} )
   sed "s/<test_num>/${tnum}/"  $srcdir/test_fms/data_override/input.nml_base > input.nml
}

teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "Test 1: Cubic-Grid" {
    skip "The input file is missing"
   cp -r $srcdir/test_fms/data_override/INPUT $builddir/test_fms/data_override/INPUT
   run mpirun -n 2 ./test_data_override
   [ "$status" -eq 0 ]
}

@test "Test 2: Latlon-Grid" {
   skip "The input file is missing"
   run mpirun -n 2 ./test_data_override
   [ "$status" -eq 0 ]
   chmod 755 $builddir/test_fms/data_override/INPUT
   rm -rf $builddir/test_fms/data_override/INPUT
}
