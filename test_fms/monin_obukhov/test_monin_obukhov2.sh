#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/monin_obukhov directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test.
cp $top_srcdir/test_fms/monin_obukhov/input_base.nml input.nml

# Run test.
mpirun -n 2 ./test_monin_obukhov
