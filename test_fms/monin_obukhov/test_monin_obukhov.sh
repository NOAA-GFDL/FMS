#!/bin/sh

export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/monin_obukhov/input.nml input.nml

bats $srcdir/test_fms/monin_obukhov/test_monin_obukhov.bats
