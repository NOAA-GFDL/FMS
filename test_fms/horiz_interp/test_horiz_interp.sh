#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/horiz_interp/input.nml input.nml

bats $srcdir/test_fms/horiz_interp/test_horiz_interp.bats
