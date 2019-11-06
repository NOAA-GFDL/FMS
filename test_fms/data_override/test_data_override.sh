#!/bin/sh

# Cause test to fail if any command fails.
set -e

cp $srcdir/test_fms/data_override/data_table_base data_table
cp $srcdir/test_fms/data_override/diag_table_base diag_table

bats $srcdir/test_fms/data_override/test_data_override.bats
