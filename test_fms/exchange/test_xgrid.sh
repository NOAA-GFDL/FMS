#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/exchange/input.nml input.nml

bats $srcdir/test_fms/exchange/test_xgrid.bats
