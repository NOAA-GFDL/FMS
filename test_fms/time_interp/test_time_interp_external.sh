#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

bats $srcdir/test_fms/time_interp/test_time_interp_external.bats
