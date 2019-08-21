#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/drifters/input.nml input.nml

bats $srcdir/test_fms/drifters/test_drifters_io.bats
