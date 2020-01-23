#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Setup the run directory
##tnum=$( printf "%2.2d" ${BATS_TEST_NUMBER} )
##rm -f diag_test_${tnum}* > /dev/null 2>&1
##sed "s/<test_num>/${tnum}/" input.nml_base > input.nml

run_test test_mpp_pset 2 skip
