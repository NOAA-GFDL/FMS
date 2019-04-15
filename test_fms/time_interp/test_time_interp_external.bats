@test "1" {
    skip "This test is going to fail"
    run mpirun -n 6 ./test_time_interp_ext
    [ "$status" -eq 0 ]
}
