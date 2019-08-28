#!/bin/sh

cp $srcdir/test_fms/mpp/input.nml_base input.nml_base

bats $srcdir/test_fms/mpp/test_mpp_domains.bats
