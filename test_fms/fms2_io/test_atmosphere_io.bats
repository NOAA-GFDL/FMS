setup() {
   # Clear out any previous test files
   rm -f *.nc > /dev/null 2>&1
   rm -f *.nc.* > /dev/null 2>&1
   rm -f logfile.*.out > /dev/null 2>&1
   rm -f input.nml
   ln -s ${srcdir}/test_fms/fms2_io/input.nml input.nml
}

@test "1" {
   skip
   run mpirun -n 6 ./test_atmosphere_io
   [ "$status" -eq 0 ]
}
