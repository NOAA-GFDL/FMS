#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/interpolator/input.nml input.nml
cp $srcdir/test_fms/interpolator/diag_table diag_table

bats $srcdir/test_fms/interpolator/test_interpolator.bats
