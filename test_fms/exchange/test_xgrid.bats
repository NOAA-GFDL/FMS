@test "1" {
    skip "The input file is missing"
   cp -r $srcdir/test_fms/exchange/INPUT INPUT

   run mpirun -n 36 ./test_xgrid
   [ "$status" -eq 0 ]

   chmod 755 $builddir/test_fms/exchange/INPUT 
   rm -rf $builddir/test_fms/exchange/INPUT
}
