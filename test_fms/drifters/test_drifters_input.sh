#!/bin/sh

ncgen -o test_drifters_input.nc $srcdir/test_fms/drifters/drifters_inp_test_3d.cdl
bats $srcdir/test_fms/drifters/test_drifters_input.bats
