#!/bin/sh

# Cause test to fail if any command fails.
set -e

cp $srcdir/test_fms/axis_utils/input_base.nml input.nml

mpirun -n 2 ./test_axis_utils
