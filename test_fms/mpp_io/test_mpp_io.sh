#!/bin/sh

cp $srcdir/test_fms/mpp_io/input_base.nml input.nml

bats $srcdir/test_fms/mpp_io/test_mpp_io.bats
