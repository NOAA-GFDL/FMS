#!/bin/sh

# Cause test to fail if any command fails.
set -e

cp $srcdir/test_fms/axis_utils/input_base.nml input.nml

bats $srcdir/test_fms/axis_utils/test_axis_utils.bats
