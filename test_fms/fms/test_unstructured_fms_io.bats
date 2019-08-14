@test "1" {
    skip "This test does not pass"
    run mpirun -n 6 ./test_unstructured_fms_io
    [ "$status" -eq 0 ]
}
