@test "1" {
    run mpirun -n 6 ./test_mpi_domains
    [ "$status" -eq 0 ]
}
