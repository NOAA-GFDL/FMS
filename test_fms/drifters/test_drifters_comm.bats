@test "1" {
    run mpirun -n 6 ./test_drifters_comm
    [ "$status" -eq 0 ]
}
