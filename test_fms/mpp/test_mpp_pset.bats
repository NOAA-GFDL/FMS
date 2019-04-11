@test "1" {
    run mpirun -n 6 ./test_mpp_pset
    [ "$status" -eq 0 ]
}
