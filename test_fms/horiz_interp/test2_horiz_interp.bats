@test "1" {
    skip "The input files are missing"
    run mpirun -n 6 ./test2_horiz_interp
    [ "$status" -eq 0 ]
}
