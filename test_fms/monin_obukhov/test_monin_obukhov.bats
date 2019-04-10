@test "1" {
    run mpirun -n 6 ./test_monin_obukhov
    [ "$status" -eq 0 ]
}
