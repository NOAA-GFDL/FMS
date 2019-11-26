#!/bin/sh

# Cause test to fail if any command fails.
set -e

bats $srcdir/test_fms/diag_manager/test_diag_manager.bats
