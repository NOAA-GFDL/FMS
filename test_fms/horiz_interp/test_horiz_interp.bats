@test "1" {
    run mpirun -n 6 ./test_horiz_interp
    [ "$status" -eq 0 ]
}
