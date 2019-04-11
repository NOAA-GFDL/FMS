@test "1" {
    run mpirun -n 6 ./test_time_interp
    [ "$status" -eq 0 ]
}
