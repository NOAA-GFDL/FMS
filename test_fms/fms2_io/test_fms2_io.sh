#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms2_io directory.

# Authors: Raymond Menzel
# Jessica Liptak

# Set common test settings.
. ../test_common.sh
# create a new input namelist from template
cp ${srcdir}/test_fms/fms2_io/input_base.nml input.nml
# run the tests
run_test test_fms2_io 6 skip
