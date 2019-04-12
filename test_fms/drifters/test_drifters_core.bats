@ test "1" {
    run mpirun -n 6 ./test_drifters_core
    [ "$status" -eq 0 ]
}
