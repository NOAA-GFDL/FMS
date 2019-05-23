@test "1" {
   run mpirun -n 32 ./test_mpp_io
   [ "$status" -eq 0 ]
}
