#!/bin/sh

cp $srcdir/test_fms/exchange/input_base.nml input.nml

bats $srcdir/test_fms/exchange/test_xgrid.bats
