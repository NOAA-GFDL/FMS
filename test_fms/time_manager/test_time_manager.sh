#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/time_manager/input.nml input.nml

bats $srcdir/test_fms/time_manager/test_time_manager.bats
