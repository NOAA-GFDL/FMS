#!/bin/sh
export PATH="$PATH:../bats/bin"
bats test_time_interp.bats
