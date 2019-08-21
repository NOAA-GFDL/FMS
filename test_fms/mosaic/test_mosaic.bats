@test "1" {
    skip "The input files are missing"
    run mpirun -n 6 ./test_mosaic
    [ "$status" -eq 0 ]
}
