#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp $srcdir/test_fms/mpp/input.nml_base input.nml_base

bats $srcdir/test_fms/mpp/test_mpp_domains.bats
