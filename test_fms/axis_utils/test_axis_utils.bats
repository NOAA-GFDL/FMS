@test "1" {
    run mpirun -n 6 ./test_axis_utils
    [ "$status" -eq 0 ]
}
