#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test.
cp $top_srcdir/test_fms/mpp_io/input_base.nml input.nml

if test "$mpi_launcher" != "" ; then 
    npes="-n 1"
fi

$mpi_launcher $npes ./test_mpp_io
