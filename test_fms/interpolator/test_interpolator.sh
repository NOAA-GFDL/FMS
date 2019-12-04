#!/bin/sh

cp $srcdir/test_fms/interpolator/input_base.nml input.nml
cp $srcdir/test_fms/interpolator/diag_table_base diag_table

bats $srcdir/test_fms/interpolator/test_interpolator.bats
