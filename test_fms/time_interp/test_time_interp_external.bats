@test "1" {
    run mpirun -n 6 ./test_time_interp_external
    [ "$status" -eq 0 ]
}
