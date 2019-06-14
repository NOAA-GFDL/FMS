@test "1" {
   run mpirun -n 8 ./test_mpp
   [ "$status" -eq 0 ]
}
