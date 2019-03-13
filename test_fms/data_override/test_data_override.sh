#!/bin/sh
export PATH="$PATH:../bats/bin"
bats test_data_override.bats
