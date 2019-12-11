#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/field_manager directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# If there is a mpi launcher set the number of processors to 2, otherwise just ./
if test "$mpi_launcher" != "" ; then 
    npes="-n 2"
fi

# Copy files for test.
cp $top_srcdir/test_fms/interpolator/input_base.nml input.nml
cp $top_srcdir/test_fms/interpolator/diag_table_base diag_table

# Test is skipped in bats file.
#$mpi_launcher $npes ./test_interpolator
