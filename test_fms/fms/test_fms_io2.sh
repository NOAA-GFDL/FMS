#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# These tests are skipped in bats files.
run_test test_fms_io 2 skip
run_test test_unstructured_fms_io 2 skip
