#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

bats $srcdir/test_fms/drifters/test_drifters_core.bats
