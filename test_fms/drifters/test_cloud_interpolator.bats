@test "1" {
    run mpirun -n 6 ./test_cloud_interpolator
    [ "$status" -eq 0 ]
}
