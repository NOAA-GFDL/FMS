#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/mosaic/input_nml input_nml

bats $srcdir/test_fms/mosaic/test_mosaic.bats
