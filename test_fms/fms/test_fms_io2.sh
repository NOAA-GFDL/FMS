#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy files for test.

# These tests are skipped in bats files.
#mpirun -n 2 ./test_fms_io
#mpirun -n 2 ./test_unstructured_fms_io
