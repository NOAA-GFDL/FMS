#!/bin/sh
# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms2_io directory.

# Set common test settings.
. ../test_common.sh

# link to the input namelist
cp ${srcdir}/test_fms/fms2_io/input_base.nml input.nml

# run the tests
run_test test_atmosphere_io 6 skip
