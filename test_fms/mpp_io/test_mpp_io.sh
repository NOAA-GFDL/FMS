#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/mpp_io/input.nml input.nml

bats $srcdir/test_fms/mpp_io/test_mpp_io.bats
