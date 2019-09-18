teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
    run mpirun -n 1 ./test_time_manager
    [ "$status" -eq 0 ]
}
