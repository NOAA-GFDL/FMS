@test "1" {
   run mpirun -n 36 ./test_xgrid
   [ "$status" -eq 0 ]
}
