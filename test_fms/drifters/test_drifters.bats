@test "1" {
    run mpirun -n 6 ./test_drifters
    [ "$status" -eq 0 ]
}
