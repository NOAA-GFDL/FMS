#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/axis_utils/input.nml input.nml

bats $srcdir/test_fms/axis_utils/test_axis_utils.bats
