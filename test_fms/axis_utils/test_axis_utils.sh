#!/bin/sh

cp $srcdir/test_fms/axis_utils/input.nml input.nml

bats $srcdir/test_fms/axis_utils/test_axis_utils.bats
