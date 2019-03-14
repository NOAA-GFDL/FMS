@test "1" {
   run mpirun -n 6 ./test_data_override
   [ "$status" -eq 0 ]
}
