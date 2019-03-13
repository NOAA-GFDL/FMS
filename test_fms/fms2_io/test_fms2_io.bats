setup() {
   # Clear out any previous test files
   rm -f *.nc > /dev/null 2>&1
   rm -f *.nc.* > /dev/null 2>&1
   rm -f logfile.*.out > /dev/null 2>&1
}

@test "1" {
   run mpirun -n 6 ./test_fms2_io
   [ "$status" -eq 0 ]
}
