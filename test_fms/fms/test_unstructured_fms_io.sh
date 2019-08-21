#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

bats $srcdir/test_fms/fms/test_unstructured_fms_io.bats
