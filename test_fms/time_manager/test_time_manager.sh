#!/bin/sh

cp $srcdir/test_fms/time_manager/input.nml input.nml

bats $srcdir/test_fms/time_manager/test_time_manager.bats
