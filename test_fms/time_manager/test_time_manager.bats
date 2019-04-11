@test "1" {
    run mpirun -n 6 ./test_time_manager
    [ "$status" -eq 0 ]
}
