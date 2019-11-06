teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
    skip "The input file is missing"
    run ./test_drifters
    [ "$status" -eq 0 ]
}
