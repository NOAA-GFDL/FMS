#!/bin/sh

cp $srcdir/test_fms/horiz_interp/input.nml input.nml

bats $srcdir/test_fms/horiz_interp/test_horiz_interp.bats
