#!/bin/sh
#
# author @underwoo
#
#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
#***********************************************************************
#
# Script to test the mpp_memutils_mod Fortran module code.

# Set common test settings
. ../test-lib.sh

# All tests use blank input.nml file.
touch input.nml

test_expect_success "begin/end routines of mpp_memutils_mod" '
    mpirun -n 1 ./test_mpp_memutils_begin_end
'

test_expect_success "mpp_print_memuse_stats" '
    mpirun -n 1 ./test_mpp_print_memuse_stats_stderr
'

test_expect_success "mpp_print_memuse_stats to file (stdout)" '
    mpirun -n 1 ./test_mpp_print_memuse_stats_file
'

test_expect_failure "failure caught if mpp_memuse_begin called multiple times" '
    mpirun -n 1 ./test_mpp_memutils_begin_2x
'

test_expect_failure "failure caught if mpp_memuse_end called before mpp_memuse_begin" '
    mpirun -n 1 ./test_mpp_memutils_end_before_begin
'

test_done
