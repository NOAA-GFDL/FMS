#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

bats $srcdir/test_fms/horiz_interp/test2_horiz_interp.bats
