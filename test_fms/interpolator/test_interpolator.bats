@test "1" {
    run mpirun -n 6 ./test_interpolator
    [ "$status" -eq 0 ]
}
