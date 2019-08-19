#!/bin/sh
export PATH="$PATH:$srcdir/test_fms/bats/bin"
export PATH="$PATH:$srcdir/test_fms/bats/libexec"

cp -r $srcdir/test_fms/data_override/INPUT INPUT
cp $srcdir/test_fms/data_override/data_table data_table
cp $srcdir/test_fms/data_override/diag_table diag_table

bats $srcdir/test_fms/data_override/test_data_override.bats
