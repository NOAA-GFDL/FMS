@test "1" {
    run mpirun -n 6 ./test_quicksort
    [ "$status" -eq 0 ]
}
