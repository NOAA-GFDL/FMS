#!/bin/sh

cp $srcdir/test_fms/field_manager/input_base.nml input.nml
cp $srcdir/test_fms/field_manager/field_table_base field_table

bats $srcdir/test_fms/field_manager/test_field_manager.bats
