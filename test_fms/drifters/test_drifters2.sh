#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/drifters directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test_drifters_io test.
cp $top_srcdir/test_fms/drifters/input_base.nml input.nml

# Run tests.

echo "1: Test_drifters_io"
run_test test_drifters_io 2

echo "2: Test_cloud_interpolator"
run_test test_cloud_interpolator 2

echo "3: Test_drifters_comm"
run_test test_drifters_comm 1

echo "4: Test_drifters_core"
run_test test_drifters_core 2

echo "5: Test_quicksort"
run_test test_quicksort 2
