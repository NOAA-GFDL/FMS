#!/bin/sh

cp $srcdir/test_fms/fft/input.nml input.nml

bats $srcdir/test_fms/fft/test_fft.bats
