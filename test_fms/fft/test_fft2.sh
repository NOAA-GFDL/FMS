#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fft directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test.
cp $top_srcdir/test_fms/fft/input_base.nml input.nml

# Test is skipped in bats file.
#mpirun -n 2 ./test_fft
