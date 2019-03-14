@test "1" {
   run mpirun -n 6 ./test_mpp
   [ "$status" -eq 0 ]
}
