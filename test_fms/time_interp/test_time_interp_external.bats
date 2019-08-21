@test "1" {
    skip "The input files are missing"
    run mpirun -n 6 ./test_time_interp_external
    [ "$status" -eq 0 ]
}
