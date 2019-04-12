@test "1" {
    run mpirun -n 6 ./test_fms_io
    [ "$status" -eq 0 ]
}
