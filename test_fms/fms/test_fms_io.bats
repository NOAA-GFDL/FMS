@test "1" {
    skip "Missing restart files"
    run mpirun -n 6 ./test_fms_io
    [ "$status" -eq 0 ]
}
