#!/bin/sh

cp $srcdir/test_fms/drifters/input_base.nml input.nml

bats $srcdir/test_fms/drifters/test_drifters_io.bats
