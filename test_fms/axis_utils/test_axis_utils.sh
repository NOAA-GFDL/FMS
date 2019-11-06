#!/bin/sh

cp $srcdir/test_fms/axis_utils/input_base.nml input.nml

bats $srcdir/test_fms/axis_utils/test_axis_utils.bats
