#!/bin/sh
export PATH="$PATH:../bats/bin"
bats test_drifters_io.bats
