@test "1" {
    run mpirun -n 6 ./test_fft
    [ "$status" -eq 0 ]
}
