@test "1" {
    run mpirun -n 6 ./test_field_manager
    [ "$status" -eq 0 ]
}
