#!/bin/sh
# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms2_io directory.

# Set common test settings.
. ../test_common.sh

# run the tests
run_test test_atmosphere_io 6 skip
