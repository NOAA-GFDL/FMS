#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/time_interp/input.nml input.nml

bats $srcdir/test_fms/time_interp/test_time_interp.bats
