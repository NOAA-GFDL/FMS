setup () {
  if [ "x$(uname -s)" = "xDarwin" ] || [ "x$TRAVIS" = "xtrue" ]
  then
    skip_test=true
  else
    skip_test=false
  fi
}

teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "MPP_IO runs with single MPI processes" {
   run mpirun -n 1 ./test_mpp_io
   [ "$status" -eq 0 ]
}

@test "MPP_IO runs with multiple MPI processes" {
    if [ "$skip_test" = "true" ]
    then
      skip "Test not reliable on Darwin"
    fi
   run mpirun -n 2 ./test_mpp_io
   [ "$status" -eq 0 ]
}
