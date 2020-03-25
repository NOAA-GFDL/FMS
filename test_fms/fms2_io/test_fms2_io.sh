#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms2_io directory.

# Authors: Raymond Menzel
# Jessica Liptak

# Set common test settings.
. ../test_common.sh
# Clear out any previous test files
rm -f *.nc > /dev/null 2>&1
rm -f *.nc.* > /dev/null 2>&1
rm -f logfile.*.out > /dev/null 2>&1
# run the tests
run_test test_fms2_io 6
