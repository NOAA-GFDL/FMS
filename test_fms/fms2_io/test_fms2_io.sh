#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms2_io directory.

# Authors: Raymond Menzel
# Jessica Liptak
#
# Set common test settings.
. ../test_common.sh
# make a dummy file for mpp_init to read
printf "EOF\n&dummy\nEOF" | cat > input.nml
# run the tests
run_test test_fms2_io 6 skip
