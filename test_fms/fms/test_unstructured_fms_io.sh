#!/bin/sh
export PATH="$PATH:../bats/bin"
bats test_unstructured_fms_io.bats
