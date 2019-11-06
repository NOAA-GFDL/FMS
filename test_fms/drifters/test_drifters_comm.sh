#!/bin/sh

# Cause test to fail if any command fails.
set -e

bats $srcdir/test_fms/drifters/test_drifters_comm.bats
