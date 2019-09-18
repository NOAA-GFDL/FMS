teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
    run mpirun -n 2 ./test_field_manager
    [ "$status" -eq 0 ]
}
