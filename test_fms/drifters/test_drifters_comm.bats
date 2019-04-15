@test "1" {
    skip "This test is going to fail"
    run mpirun -n 6 ./test_drifters_comm
    [ "$status" -eq 0 ]
}
