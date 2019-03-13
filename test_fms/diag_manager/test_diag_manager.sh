#!/bin/sh
export PATH="$PATH:../bats/bin"
bats test_diag_manager.bats
