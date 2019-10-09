#!/bin/sh

cp $srcdir/test_fms/mpp/input_base.nml input.nml_base

bats $srcdir/test_fms/mpp/test_mpp_domains.bats
