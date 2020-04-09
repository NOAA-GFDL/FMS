#!/bin/sh
# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms2_io directory.

# Set common test settings.
. ../test_common.sh
# make an input.nml for mpp_init to read
printf "EOF\n&dummy\nEOF" | cat > input.nml
# run the tests
run_test test_atmosphere_io 6 skip
