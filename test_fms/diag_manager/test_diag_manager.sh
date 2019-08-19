#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"
bats $srcdir/test_fms/diag_manager/test_diag_manager.bats
