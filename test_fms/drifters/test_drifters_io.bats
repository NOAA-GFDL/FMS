@test  "1" {
    run mpirun -n 6 ./test_drifters_io
    [ "$status" -eq 0 ]
}
