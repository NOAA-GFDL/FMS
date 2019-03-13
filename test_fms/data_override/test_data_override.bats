@test "1" {
   skip "TODO: get this test to pass"
   run mpirun -n 6 ./test_data_override
   [ "$status" -eq 0 ]
}
