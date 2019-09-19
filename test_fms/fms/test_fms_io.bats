teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
    skip "Missing restart files"
    run mpirun -n 2 ./test_fms_io
    [ "$status" -eq 0 ]
}
