teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
   run  ./test_drifters_comm
    [ "$status" -eq 0 ]
}
