@test "1" {
    run mpirun -n 6 ./test_drifters_input
    [ "$status" -eq 0 ]
}
