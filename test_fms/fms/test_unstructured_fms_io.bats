@ test "1" {
    run mpirun -n 6 ./test_unstructured_fms_io
    [ "$status" -eq 0 ]
}
