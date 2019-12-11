#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/time_manager directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test.
cp $top_srcdir/test_fms/time_manager/input_base.nml input.nml

# If there is a mpi launcher set the number of processors to 2, otherwise just ./
if test "$mpi_launcher" != "" ; then 
    npes="-n 1"
fi

$mpi_launcher $npes ./test_time_manager
