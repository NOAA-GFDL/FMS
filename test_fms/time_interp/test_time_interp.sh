#!/bin/sh

cp $srcdir/test_fms/time_interp/input_base.nml input.nml

bats $srcdir/test_fms/time_interp/test_time_interp.bats
