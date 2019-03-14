@test "1" {
   skip "TODO: get this to pass"
   run mpirun -n 1 ./test_mpp_io
   [ "$status" -eq 0 ]
}
