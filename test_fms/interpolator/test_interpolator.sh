#!/bin/sh

cp $srcdir/test_fms/interpolator/input.nml input.nml
cp $srcdir/test_fms/interpolator/diag_table diag_table

bats $srcdir/test_fms/interpolator/test_interpolator.bats
