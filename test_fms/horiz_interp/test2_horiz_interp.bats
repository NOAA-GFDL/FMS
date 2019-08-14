@test "1" {
    run mpirun -n 6 ./test2_horiz_interp
    [ "$status" -eq 0 ]
}
