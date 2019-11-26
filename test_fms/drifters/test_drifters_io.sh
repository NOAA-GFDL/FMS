#!/bin/sh

# Cause test to fail if any command fails.
set -e

cp $srcdir/test_fms/drifters/input_base.nml input.nml

bats $srcdir/test_fms/drifters/test_drifters_io.bats
