#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/horiz_interp directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test.
cp $top_srcdir/test_fms/horiz_interp/input_base.nml input.nml

mpirun -n 2 ./test_horiz_interp

# Other test is skipped in bats.
#mpirun -n 2 ./test_horiz_interp
