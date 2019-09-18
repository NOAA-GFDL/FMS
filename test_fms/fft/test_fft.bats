teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
    skip "Fails in 32bit/Mixed mode"
    run mpirun -n 2 ./test_fft
    [ "$status" -eq 0 ]
}
