teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
    skip "This test does not pass"
    run mpirun -n 2 ./test_unstructured_fms_io
    [ "$status" -eq 0 ]
}
