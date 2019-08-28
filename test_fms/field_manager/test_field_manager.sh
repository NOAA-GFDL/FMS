#!/bin/sh

cp $srcdir/test_fms/field_manager/input.nml input.nml
cp $srcdir/test_fms/field_manager/field_table field_table

bats $srcdir/test_fms/field_manager/test_field_manager.bats
