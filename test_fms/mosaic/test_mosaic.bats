@test "1" {
    run mpirun -n 6 ./test_mosaic
    [ "$status" -eq 0 ]
}
