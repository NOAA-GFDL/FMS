#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/axis_utils directory.

# Ed Hartnett 11/26/19

# Set common test settings.
. ../test_common.sh

# Output for TAP.
echo 1..1

# Copy and rename namelist file.
cp $top_srcdir/test_fms/axis_utils/input_base.nml input.nml

# Run the test.
mpirun -n 2 ./test_axis_utils

# Output for TAP.
echo ok 1
