#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/field_manager directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy files for test.
#cp $top_srcdir/test_fms/mosaic/input_nml input_nml
run_test test_mosaic 2 skip
