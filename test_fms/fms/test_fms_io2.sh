#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

if test "$mpi_launcher" != "" ; then 
    npes="-n 2"
fi

# Copy files for test.

# These tests are skipped in bats files.
#$mpi_launcher $npes ./test_fms_io
#$mpi_launcher $npes ./test_unstructured_fms_io
