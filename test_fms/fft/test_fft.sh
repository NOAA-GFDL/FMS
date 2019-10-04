#!/bin/sh

cp $srcdir/test_fms/fft/input_base.nml input.nml

bats $srcdir/test_fms/fft/test_fft.bats
