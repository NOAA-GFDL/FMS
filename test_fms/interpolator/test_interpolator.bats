@test "1" {
    skip "The input files are missing"
    run mpirun -n 6 ./test_interpolator
    [ "$status" -eq 0 ]
}
