#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/drifters directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test_drifters_io test.
cp $top_srcdir/test_fms/drifters/input_base.nml input.nml

# Run tests.
mpirun -n 2 ./test_drifters_io
mpirun -n 2 ./test_cloud_interpolator
./test_drifters_comm
mpirun -n 2 ./test_drifters_core
mpirun -n 2 ./test_quicksort
