#!/bin/sh

cp $srcdir/test_fms/drifters/input.nml input.nml

bats $srcdir/test_fms/drifters/test_drifters_io.bats
