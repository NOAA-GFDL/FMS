#!/bin/sh
export PATH="$PATH:../bats/bin"
bats test_interpolator.bats
